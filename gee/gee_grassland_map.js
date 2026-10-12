// =============================================================================
// Grassland condition anomaly — WALL-TO-WALL MAP for a small AOI
//
// Companion to gee_grassland_extract.js, and deliberately separate:
//   * the extract script SAMPLES parcels across the canton -> statistics
//   * this script covers EVERY grassland parcel in a small box -> a map
// A scattered sample cannot be mapped; a wall-to-wall canton cannot be afforded.
//
// This is cheap. It needs two composites, not 60 dates x 8 years, so the cost
// scales with AOI area only. Keep the box small and this runs interactively.
// =============================================================================

// ----------------------------------------------------------------------------
// 1. Config
// ----------------------------------------------------------------------------
// >>> drought gee render: generated from config.yaml. Do not edit by hand —
// change config.yaml and run `drought gee render`.
// Paste-ready copy: drought --set gee.project=<project> gee render --out output/gee
var PROJECT = "YOUR_EE_PROJECT";   // Earth Engine project with the assets
var ASSETS = 'projects/' + PROJECT + '/assets/';
var MASK_ASSET = ASSETS + "grassland_mask_sg";   // uploaded from build_mask()'s GeoTIFF
var PARCEL_PREFIX = ASSETS + "nutzungsflaechen_";   // + year
var CONFIG_HASH = "bdedabe82ea4";   // settings + mask this copy was rendered with; goes into export names
var PURITY_MIN = 0.9;   // mask.purity_min
var CLEAR_THRESHOLD = 0.6;   // Cloud Score+ cs_cdf, gee.clear_threshold
var MASK_BANDS = ["mown", "grazed", "intensive", "extensive"];   // band order of the mask raster, as written by build_mask()
var MOWN_CODES = [611, 612, 613, 634, 635, 650, 694, 621, 622, 623, 935];   // landuse.mown
var GRAZED_CODES = [616, 617, 618, 660, 693, 625, 930, 933];   // landuse.grazed
var YEAR_DRY = 2026;   // project.current_year
var BASELINE_YEARS = [2019, 2020, 2021, 2022, 2023, 2024, 2025];   // project.baseline_years
// <<< end of generated block

var ACTIVE_SITE = 'rheintal';
// 'rheintal' | 'linthebene' | 'alpstein' | 'diemtigtal'

// NOTE ON THE MEDIA SITES: the 2026 water-flight coverage was Fribourg
// (Schwarzsee), Bern (Diemtigen), Schwyz and Graubuenden -- NOT St. Gallen.
// Your Nutzungsflaechen are cantonal, so an out-of-canton AOI has no parcel
// data. Sites flagged parcels:false fall back to a pixel-based grassland mask
// (ESA WorldCover) instead. Weaker, but it works anywhere.
//
// COORDINATES BELOW ARE APPROXIMATE -- check each box on the map before use.
var SITES = {
  rheintal: {
    name: 'St. Galler Rheintal - valley-floor grassland',
    geom: ee.Geometry.Rectangle([9.50, 47.32, 9.72, 47.50]),
    zoom: 12, stratum: 'lowland', parcels: true
  },
  linthebene: {
    name: 'Linthebene - intensive grassland on drained peat',
    geom: ee.Geometry.Rectangle([8.90, 47.16, 9.12, 47.30]),
    zoom: 12, stratum: 'lowland', parcels: true
  },
  alpstein: {
    name: 'Alpstein / Churfirsten - SG summering area',
    geom: ee.Geometry.Rectangle([9.25, 47.15, 9.48, 47.30]),
    zoom: 11, stratum: 'summering', parcels: true
  },
  diemtigtal: {
    name: 'Diemtigtal BE - army water flights, Aug 2026 (no parcel data)',
    geom: ee.Geometry.Rectangle([7.48, 46.52, 7.68, 46.68]),
    zoom: 11, stratum: 'summering', parcels: false
  }
};

var site = SITES[ACTIVE_SITE];
var AOI = site.geom;
print('Active site:', site.name);

var assetForYear = function(year) {
  return PARCEL_PREFIX + year;
};

// The persistence raster built by drought.mask.build_mask(). Parcel polygons come
// from the vector asset (this map is per-parcel), but the anomaly is computed
// ONLY over pixels that were grassland in every year — otherwise a parcel that
// was arable during part of the baseline gets a spurious anomaly.
// MASK_ASSET: rendered from config.yaml.
// Four float bands of PSF-weighted purity 0-1, on the 10 m Sentinel grid.
// Persistence is NOT in the raster — it is in the parcel selection that built
// the mask (grassland in every year AND the same management category), so the
// bands carry only purity and the sensor halo is intact.
// MASK_BANDS: rendered (the band order build_mask() writes).

// Threshold here, not in the mask: sweeping this costs one extraction, not a
// rebuild and re-upload. On the previous 20 m mask the NDMI anomaly moved by
// under 1% of its own size across 0.05-0.95; not yet re-run at 10 m.
// mown/grazed is about whether a CUT SIGNAL exists; intensive/extensive is
// about the FARMING SYSTEM, and the two do not nest: 'mown' pools 613 (4-6 cuts
// from early May) with 611/612 (~2 cuts, none before 15 June).
// The mask bands are CONTINUOUS purity 0-1: the fraction of each 20 m cell
// that lies inside a parcel. Threshold here — sweeping this costs one
// extraction, not a mask rebuild and re-upload.
// PURITY_MIN: rendered from config.yaml (mask.purity_min).
// NOTE on reading a float mask at a coarser scale: GEE's default pyramiding for
// float is 'mean', so at 50 m you would get the MEAN of four 20 m purities — a
// plausible-looking number that is the purity of nothing. Upload the asset with
// pyramidingPolicy {'.default': 'sample'}, or reduce at the native 20 m.
var USE_PERSISTENCE = true;   // false = every pixel inside the parcel
var GEOM_YEAR  = 2025;   // geometries + declared use taken from this year
var FIELD_CODE = 'lnf_code';

// Verified against the official LNF code list. TWO lists, not four: the split
// that matters is mown vs grazed, because grazing produces no cut signal and
// has different interannual variability. Whether a code is an LN code or a
// Soemmerung code is not a distinction this analysis acts on, and an earlier
// version that switched code lists by site was really just guessing which
// codes occur where. Filter on all of them and let the data decide -- a
// lowland box simply contains no Soemmerung parcels, and vice versa.
//
// mown: 611/612/613 Dauerwiesen, 634/635 Uferwiesen, 650 aggregate,
//       694 BFF Gruenflaeche, 621/622/623/935 Heuwiesen im Soemmerungsgebiet
// MOWN_CODES: rendered from config.yaml (landuse.mown).
// grazed: 616/617 Weiden, 618 Waldweiden, 660 aggregate, 693 BFF Weiden,
//         625 Waldweiden Soemmerung, 930 Soemmerungsweiden, 933 Gemeinschaftsweiden
// GRAZED_CODES: rendered from config.yaml (landuse.grazed).
// Deliberately EXCLUDED everywhere: 601/602 Kunstwiesen (temporary ley, arable
// rotation), 631/632 seed production, 851/936 Streueflaechen, 697/698 (regime
// cannot be determined).

// YEAR_DRY, BASELINE_YEARS: rendered from config.yaml (project).

// Window depends on the stratum. This is the ONLY thing `stratum` still
// controls, and it is a physical argument rather than a bookkeeping one:
//  lowland   1 May - 31 Aug. Long on purpose: in a narrow window a parcel cut
//            three days before the acquisition looks catastrophic. ~4 months
//            averages over one to two cut cycles and measures condition.
//  summering 15 Jul - 15 Sep. Later on purpose: at 1500-2200 m the season
//            starts at snowmelt, whose date swings by weeks between years. 2026
//            had a warm dry spring, so an early-May window would compare
//            snow-covered baseline pixels against already-green 2026 pixels and
//            produce a spurious POSITIVE anomaly. A late window sits after
//            snowmelt in both early and normal years.
var WINDOWS = {
  lowland:   {start: 121, end: 243},
  summering: {start: 196, end: 258}
};
var DOY_START = WINDOWS[site.stratum].start;
var DOY_END   = WINDOWS[site.stratum].end;

var INDEX = 'NDMI';    // 'NDMI' (canopy water, responds earlier) or 'NDVI'

var MIN_AREA_HA   = 0.25;
var EDGE_BUFFER_M = -10;

// ----------------------------------------------------------------------------
// 2. Grassland units -- parcels where declared, pixel mask otherwise
// ----------------------------------------------------------------------------
var GRASS_CODES = MOWN_CODES.concat(GRAZED_CODES);

var parcels = null;
var pixelMask = null;

if (site.parcels) {
  parcels = ee.FeatureCollection(assetForYear(GEOM_YEAR))
    .filterBounds(AOI)
    .filter(ee.Filter.inList(FIELD_CODE, GRASS_CODES))
    .map(function(f) {
      var code = ee.Number.parse(ee.Algorithms.String(f.get(FIELD_CODE)));
      return f
        .set('regime', ee.Algorithms.If(
          ee.List(MOWN_CODES).contains(code), 'mown', 'grazed'))
        .set('area_ha', f.geometry().area(1).divide(10000));
    })
    .filter(ee.Filter.gte('area_ha', MIN_AREA_HA))
    .map(function(f) { return f.setGeometry(f.geometry().buffer(EDGE_BUFFER_M, 1)); });
  print('Grassland parcels in AOI:', parcels.size());
} else {
  // Fallback: ESA WorldCover grassland (class 30), 10 m, matches the S2 grid.
  // No regime split is possible here -- mown and grazed cannot be told apart
  // without declarations, so the statistics are pooled and weaker.
  pixelMask = ee.ImageCollection('ESA/WorldCover/v200').first()
    .select('Map').eq(30).selfMask().clip(AOI);
  print('No parcel data for this site - using ESA WorldCover grassland mask.');
}

// ----------------------------------------------------------------------------
// 3. Composites
// ----------------------------------------------------------------------------
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED');
var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');

var indexComposite = function(year) {
  return s2
    .filterBounds(AOI)
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .filter(ee.Filter.calendarRange(DOY_START, DOY_END, 'day_of_year'))
    .linkCollection(csPlus, ['cs_cdf'])
    .map(function(img) {
      var m = img.updateMask(img.select('cs_cdf').gte(CLEAR_THRESHOLD));
      return INDEX === 'NDMI'
        ? m.normalizedDifference(['B8', 'B11']).rename(INDEX)
        : m.normalizedDifference(['B8', 'B4']).rename(INDEX);
    })
    .median()          // median across the season = condition, not mowing phase
    .clip(AOI);
};

var dry = indexComposite(YEAR_DRY);
var baseCol = ee.ImageCollection(BASELINE_YEARS.map(indexComposite));
var baseMed = baseCol.median();
var baseSd = baseCol.reduce(ee.Reducer.stdDev()).rename('sd');

var anomaly = dry.subtract(baseMed).rename('d');
// z-score is the better map: it says "unusual relative to this parcel's own
// year-to-year variability", which a raw difference cannot.
var z = anomaly.divide(baseSd.max(0.01)).rename('z');

// ----------------------------------------------------------------------------
// 4. Aggregation and rendering
// ----------------------------------------------------------------------------
var divPalette = ['#67001f', '#b2182b', '#ef8a62', '#fddbc7',
                  '#f7f7f7', '#d1e5f0', '#67a9cf', '#2166ac'];

// SATURATION FIX. The original +/-2.5 range clipped almost every parcel in 2026
// (Linthebene mean z was -2.4), so the map rendered uniformly dark red and all
// spatial structure was lost. Two changes:
//   * a wider z range, so extreme parcels are still distinguishable;
//   * a RAW DIFFERENCE layer, which has no small-denominator pathology at all.
// Raw difference is also the honest cross-check on the grazed-vs-mown gap: z
// divides by interannual SD, and mown parcels vary more between years because
// cut timing shifts, so a bigger denominator shrinks their z for the same real
// change. If the regime gap survives in dVis it is real; if not, it was an
// artefact of the denominator.
var Z_RANGE = 5;      // was 2.5
var D_RANGE = 0.20;   // raw index units
var zVis = {min: -Z_RANGE, max: Z_RANGE, palette: divPalette};
var dVis = {min: -D_RANGE, max: D_RANGE, palette: divPalette};

Map.centerObject(AOI, site.zoom);
Map.setOptions('SATELLITE');

var painted, stats = null;

if (site.parcels) {
  // Restrict to persistently-grassland pixels before aggregating per parcel.
  var persistMask = ee.Image(MASK_ASSET).rename(MASK_BANDS)
    .select(['mown', 'grazed']).reduce(ee.Reducer.max())
    .gte(PURITY_MIN);
  var zM    = USE_PERSISTENCE ? z.updateMask(persistMask) : z;
  var anomM = USE_PERSISTENCE ? anomaly.updateMask(persistMask) : anomaly;

  stats = zM.addBands(anomM).reduceRegions({
    collection: parcels,
    reducer: ee.Reducer.mean(),
    scale: 10,
    tileScale: 4
  }).filter(ee.Filter.notNull(['z']));

  // Choropleth by parcel -- matches the unit the statistics are computed on.
  painted = stats.reduceToImage({properties: ['z'], reducer: ee.Reducer.first()});
  var paintedD = stats.reduceToImage({properties: ['d'], reducer: ee.Reducer.first()});
  Map.addLayer(painted, zVis, 'Condition z-score (per parcel)', true);
  Map.addLayer(paintedD, dVis, 'Raw ' + INDEX + ' difference (per parcel)', false);
  Map.addLayer(parcels.style({color: '00000060', fillColor: '00000000', width: 1}),
               {}, 'Parcel boundaries', false);
} else {
  painted = z.updateMask(pixelMask);
  Map.addLayer(painted, zVis, 'Condition z-score (per pixel, WorldCover mask)', true);
}

Map.addLayer(z.updateMask(painted.mask()), zVis, 'Condition z-score (per pixel)', false);

// Categorical version -- easier to read on a slide than a continuous ramp.
var classed = ee.Image(0)
  .where(painted.lt(-2), 1)
  .where(painted.gte(-2).and(painted.lt(-1)), 2)
  .where(painted.gte(-1).and(painted.lt(1)), 3)
  .where(painted.gte(1), 4)
  .updateMask(painted.mask()).rename('class');
Map.addLayer(classed,
  {min: 1, max: 4, palette: ['#67001f', '#ef8a62', '#f7f7f7', '#2166ac']},
  'Class: severe / below / normal / above', false);

// ----------------------------------------------------------------------------
// 5. Summary numbers for the caption
// ----------------------------------------------------------------------------
if (site.parcels) {
  print('Mean parcel z-score:', stats.aggregate_mean('z'));
  ['mown', 'grazed'].forEach(function(reg) {
    var sub = stats.filter(ee.Filter.eq('regime', reg));
    print(reg + ' - n parcels:', sub.size());
    print(reg + ' - mean z:', sub.aggregate_mean('z'));
    print(reg + ' - mean raw d' + INDEX + ':', sub.aggregate_mean('d'));
    print(reg + ' - n below z = -1:', sub.filter(ee.Filter.lt('z', -1)).size());
  });
} else {
  print('Mean z-score over WorldCover grassland:', z.updateMask(pixelMask).reduceRegion({
    reducer: ee.Reducer.mean(), geometry: AOI, scale: 20,
    maxPixels: 1e8, bestEffort: true, tileScale: 4}));
  print('Fraction below z = -1:', z.lt(-1).updateMask(pixelMask).reduceRegion({
    reducer: ee.Reducer.mean(), geometry: AOI, scale: 20,
    maxPixels: 1e8, bestEffort: true, tileScale: 4}));
}

if (site.stratum === 'summering') {
  print('CAVEAT for summering sites: the binding constraint on an alp in a ' +
        'drought is often DRINKING WATER, not forage. Springs and cisterns ' +
        'running dry force early descent even where grass is adequate, and a ' +
        'vegetation index cannot see that at all. A modest anomaly here is ' +
        'fully compatible with a serious alpine crisis. Validate against ' +
        'Alpabzug dates, not against the map.');
}

// ----------------------------------------------------------------------------
// 7. Legend
// ----------------------------------------------------------------------------
var legend = ui.Panel({
  style: {position: 'bottom-left', padding: '8px 12px', backgroundColor: 'white'}
});

legend.add(ui.Label({
  value: site.name.split(' - ')[0] + ' - ' + YEAR_DRY,
  style: {fontWeight: 'bold', fontSize: '14px', margin: '0 0 2px 0'}
}));
legend.add(ui.Label({
  value: INDEX + ' z-score vs. ' + BASELINE_YEARS[0] + '-' +
         BASELINE_YEARS[BASELINE_YEARS.length - 1],
  style: {fontSize: '11px', color: '#555', margin: '0 0 6px 0'}
}));

// Continuous ramp for the z-score layer.
var rampRow = ui.Panel({layout: ui.Panel.Layout.flow('horizontal')});
divPalette.forEach(function(c) {
  rampRow.add(ui.Label({
    value: '', style: {backgroundColor: c, padding: '8px 9px', margin: '0'}
  }));
});
legend.add(rampRow);

var scaleRow = ui.Panel({
  layout: ui.Panel.Layout.flow('horizontal'),
  style: {margin: '2px 0 8px 0'}
});
['worse  -' + Z_RANGE, '0', '+' + Z_RANGE + '  better'].forEach(function(t, i) {
  scaleRow.add(ui.Label({
    value: t,
    style: {fontSize: '10px', margin: '0', stretch: i === 1 ? 'horizontal' : 'none',
            textAlign: i === 1 ? 'center' : 'left'}
  }));
});
legend.add(scaleRow);

// Class swatches, matching the categorical layer.
legend.add(ui.Label({
  value: 'Classes',
  style: {fontWeight: 'bold', fontSize: '11px', margin: '4px 0 2px 0'}
}));

var classItem = function(color, label) {
  return ui.Panel({
    widgets: [
      ui.Label({value: '', style: {
        backgroundColor: color, padding: '7px', margin: '0 6px 3px 0',
        border: '1px solid #999'}}),
      ui.Label({value: label, style: {fontSize: '11px', margin: '0 0 3px 0'}})
    ],
    layout: ui.Panel.Layout.flow('horizontal')
  });
};

legend.add(classItem('#67001f', 'Severe  (z < -2)'));
legend.add(classItem('#ef8a62', 'Below normal  (-2 to -1)'));
legend.add(classItem('#f7f7f7', 'Normal  (-1 to +1)'));
legend.add(classItem('#2166ac', 'Above normal  (z > +1)'));

legend.add(ui.Label({
  // Name the mapped class the same way the deck and the trajectory figure do.
  // The earlier wording ("mown and grazed shown together") described the
  // internals rather than the class, and read as if it contradicted the slide
  // text calling it permanent grassland. Both were true, but only one is useful
  // to a reader.
  value: site.parcels
    ? 'Permanent grassland (Dauerwiesen and Weiden), persistent 2017-2025.'
    : 'Pixel-based (ESA WorldCover) - no parcel data for this site.',
  style: {fontSize: '9px', color: '#777', margin: '6px 0 0 0'}
}));

Map.add(legend);

// ----------------------------------------------------------------------------
// 8. Export
// ----------------------------------------------------------------------------
Export.image.toDrive({
  image: painted.visualize(zVis),
  description: 'grassland_zscore_' + ACTIVE_SITE + '_' + YEAR_DRY,
  folder: 'GEE_Exports', crs: 'EPSG:2056', scale: 10,
  region: AOI, maxPixels: 1e9
});


// The parcel table, so the map and the statistics come from the same object.
if (site.parcels) {
  Export.table.toDrive({
    collection: stats,
    description: 'grassland_parcels_' + ACTIVE_SITE + '_' + YEAR_DRY,
    folder: 'GEE_Exports', fileFormat: 'CSV',
    selectors: ['area_ha', 'regime', FIELD_CODE, 'z', 'd']
  });
}
