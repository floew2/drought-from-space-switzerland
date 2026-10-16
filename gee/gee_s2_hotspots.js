// =============================================================================
// Drought 2026 — Sentinel-2 zoom-ins for selected hotspots
// Run ONCE PER SITE (set ACTIVE_SITE). For looking, not for data.
//
// Why Sentinel-2:
//   * 10 m — individual fields, shorelines, field boundaries.
//   * A STABLE overpass. Sentinel-2's mean local descending node is ~10:30 and
//     is actively maintained, so there is no solar-zenith drift between years.
//     Whatever difference you see is much closer to real surface change.
//   * All three sites are flat, so topographic shading is largely absent.
//
// What you give up, and it must be stated:
//   * The baseline is only 7 years (2019-2025). Sentinel-2 L2A does not reach
//     back far enough for a 20- or 30-year climatology.
//   * That short baseline CONTAINS the 2022 and 2018 Swiss drought summers, so
//     the "normal" is drier than a true long-term normal. The 2026 anomaly
//     computed here is therefore conservative — a floor, not a ceiling.
//   * S2A/S2B/S2C are different units; inter-sensor differences are small but
//     non-zero across a 7-year baseline.
// =============================================================================

// ----------------------------------------------------------------------------
// 1. Configuration
// ----------------------------------------------------------------------------
// >>> drought gee render: generated from config.yaml. Do not edit by hand —
// change config.yaml and run `drought gee render`.
// Paste-ready copy: drought --set gee.project=<project> gee render --out output/gee
var CLEAR_THRESHOLD = 0.6;   // Cloud Score+ cs_cdf, gee.clear_threshold
var YEAR_DRY = 2026;   // project.current_year
var BASELINE_YEARS = [2019, 2020, 2021, 2022, 2023, 2024, 2025];   // project.baseline_years
// <<< end of generated block

var ACTIVE_SITE = 'untersee';   // 'untersee' | 'rheintal' | 'linthebene'

// YEAR_DRY, BASELINE_YEARS: rendered from config.yaml (project).

var DOY_START = 166;   // mid-June to end of July: the onset of the drought
var DOY_END   = 212;

var QA_BAND        = 'cs_cdf';  // 'cs' is more instantaneous; 'cs_cdf' more stable
// CLEAR_THRESHOLD: rendered from config.yaml (gee.clear_threshold); catalogue guidance 0.50-0.65
var MAX_SCENE_CLOUD = 80;       // drop hopeless scenes early, saves compute

var WATER_METHOD = 'both';   // 'ndwi' | 'dynamicworld' | 'both'
// NDWI > 0 is a threshold someone chose. Dynamic World is a convolutional
// classifier trained on a large labelled dataset, run per Sentinel-2 scene at
// 10 m, giving a per-pixel WATER PROBABILITY instead of a hand-set cut-off.
// It is not a "foundation model" — it is a supervised land-cover CNN — but it
// is the closest thing GEE serves, and it handles turbid and shallow margins
// better than a ratio does. Where the two disagree is exactly the shallow
// shoreline this figure is about, so 'both' is the honest default.
var DW_THRESHOLD = 0.5;      // water probability counted as water
var DW_CHANGE_MIN = 0.10;    // hide changes smaller than this in the change
                             // layer: below it the difference is model noise,
                             // and painting it covers the basemap for nothing

var SHOW_ANALYSIS = true;       // NDVI anomaly + (for lakes) water-area change

// ----------------------------------------------------------------------------
// 2. Sites
// ----------------------------------------------------------------------------
// NOTE on Untersee: the Swiss shore here is canton THURGAU, not St. Gallen.
// Worth getting right on a slide. The Untersee is shallow and unregulated, so
// a modest drop in level exposes a disproportionately wide shoreline — which is
// exactly what makes it photogenic at 10 m.
var SITES = {
  untersee: {
    name: 'Untersee (Bodensee) — TG / DE',
    geom: ee.Geometry.Rectangle([8.83, 47.60, 9.22, 47.76]),
    zoom: 12,
    water: true
  },
  rheintal: {
    name: 'St. Galler Rheintal — agriculture',
    geom: ee.Geometry.Rectangle([9.50, 47.32, 9.72, 47.50]),
    zoom: 12,
    water: false
  },
  linthebene: {
    name: 'Linthebene — intensive agriculture on drained peat',
    geom: ee.Geometry.Rectangle([8.90, 47.16, 9.12, 47.30]),
    zoom: 12,
    water: false
  }
};

var site = SITES[ACTIVE_SITE];
var aoi = site.geom;
Map.centerObject(aoi, site.zoom);
Map.setOptions('SATELLITE');
print('Active site:', site.name);

// ----------------------------------------------------------------------------
// 3. Collections
// ----------------------------------------------------------------------------
// S2_SR_HARMONIZED, not S2_SR: the harmonized collection removes the -1000
// radiometric offset introduced with processing baseline 04.00 in Jan 2022.
// Without it, a baseline spanning 2019-2025 would contain a step change.
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED');
var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');

var RGB   = ['B4', 'B3', 'B2'];
var BANDS = ['B2', 'B3', 'B4', 'B8', 'B11', 'B12'];

var yearlyComposite = function(year) {
  var col = s2
    .filterBounds(aoi)
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .filter(ee.Filter.calendarRange(DOY_START, DOY_END, 'day_of_year'))
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', MAX_SCENE_CLOUD))
    .linkCollection(csPlus, [QA_BAND])
    .map(function(img) {
      return img.updateMask(img.select(QA_BAND).gte(CLEAR_THRESHOLD))
        .select(BANDS);
    });
  return col.median().clip(aoi).set('year', year);
};

var sceneCount = function(year) {
  return s2.filterBounds(aoi)
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .filter(ee.Filter.calendarRange(DOY_START, DOY_END, 'day_of_year'))
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', MAX_SCENE_CLOUD))
    .size();
};

var imgDry = yearlyComposite(YEAR_DRY);
var baseComposites = BASELINE_YEARS.map(yearlyComposite);
var imgBase = ee.ImageCollection(baseComposites).median().clip(aoi);

print('Scenes, ' + YEAR_DRY + ':', sceneCount(YEAR_DRY));
print('Scenes, ' + BASELINE_YEARS[0] + ':', sceneCount(BASELINE_YEARS[0]));
print('If the ' + YEAR_DRY + ' count is very low, widen DOY_START/DOY_END ' +
      'before trusting the composite.');

// ----------------------------------------------------------------------------
// 4. Visualisation — identical pinned stretch on both panels
// ----------------------------------------------------------------------------
var visTrue  = {bands: RGB, min: 0, max: 3000, gamma: 1.2};
var visFalse = {bands: ['B8', 'B4', 'B3'], min: 0, max: 4500, gamma: 1.2};

var baseLabel = BASELINE_YEARS[0] + '-' + BASELINE_YEARS[BASELINE_YEARS.length - 1] +
                ' median (' + BASELINE_YEARS.length + ' yr)';

Map.addLayer(imgBase, visTrue, 'Baseline ' + baseLabel + ' — true colour', true);
Map.addLayer(imgDry,  visTrue, YEAR_DRY + ' — true colour', true);
Map.addLayer(imgBase, visFalse, 'Baseline — false colour (NIR)', false);
Map.addLayer(imgDry,  visFalse, YEAR_DRY + ' — false colour (NIR)', false);

// ----------------------------------------------------------------------------
// 5. Analysis
// ----------------------------------------------------------------------------
if (SHOW_ANALYSIS) {
  var ndvi = function(img) { return img.normalizedDifference(['B8', 'B4']).rename('NDVI'); };
  var ndwi = function(img) { return img.normalizedDifference(['B3', 'B8']).rename('NDWI'); };

  // ESA WorldCover, 10 m — matches the Sentinel-2 grid, so no rescaling.
  var wc = ee.ImageCollection('ESA/WorldCover/v200').first().select('Map').clip(aoi);
  var isAgri = wc.eq(30).or(wc.eq(40));   // grassland + cropland

  var ndviBaseCol = ee.ImageCollection(baseComposites.map(ndvi));
  var anomaly = ndvi(imgDry).subtract(ndviBaseCol.median()).rename('dNDVI');
  var zscore = anomaly.divide(ndviBaseCol.reduce(ee.Reducer.stdDev()).max(0.01)).rename('z');

  Map.addLayer(anomaly.updateMask(isAgri),
    {min: -0.3, max: 0.3,
     palette: ['#8c2d04', '#fdae61', '#ffffbf', '#a6d96a', '#1a9850']},
    'NDVI anomaly on grassland+cropland', false);
  Map.addLayer(zscore.updateMask(isAgri),
    {min: -3, max: 3, palette: ['#67001f', '#f7f7f7', '#053061']},
    'NDVI z-score', false);

  var agriStats = anomaly.addBands(zscore.lt(-1).rename('below1sd'))
    .updateMask(isAgri).reduceRegion({
      reducer: ee.Reducer.mean(), geometry: aoi, scale: 20,
      maxPixels: 1e8, bestEffort: true, tileScale: 4
    });
  print('Agricultural land — mean dNDVI and fraction below -1 sigma:', agriStats);

  // --- Dynamic World water probability -------------------------------------
  // Same window, same compositing rule as the reflectance panels.
  var dwWater = function(year) {
    return ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1')
      .filterBounds(aoi)
      .filter(ee.Filter.calendarRange(year, year, 'year'))
      .filter(ee.Filter.calendarRange(DOY_START, DOY_END, 'day_of_year'))
      .select('water')
      .median()
      .clip(aoi);
  };

  // --- water extent, lake sites only ---------------------------------------
  if (site.water) {
    var px = ee.Image.pixelArea().divide(1e6);   // km2

    var pairs = [];
    if (WATER_METHOD === 'ndwi' || WATER_METHOD === 'both') {
      pairs.push({tag: 'NDWI>0', dry: ndwi(imgDry).gt(0),
                  base: ndwi(imgBase).gt(0)});
    }
    if (WATER_METHOD === 'dynamicworld' || WATER_METHOD === 'both') {
      var dwCur = dwWater(YEAR_DRY);
      var dwBase = ee.ImageCollection(BASELINE_YEARS.map(dwWater)).median();
      pairs.push({tag: 'DynamicWorld p>' + DW_THRESHOLD,
                  dry: dwCur.gte(DW_THRESHOLD),
                  base: dwBase.gte(DW_THRESHOLD)});
      Map.addLayer(dwCur, {min: 0, max: 1,
        palette: ['#000000', '#2166ac']}, YEAR_DRY + ' - DW water probability',
        false);
      // Mask the near-zero middle so unchanged ground is transparent rather
      // than a white sheet over the basemap. Only pixels that actually changed
      // are painted, which is what the layer is for.
      var dwDelta = dwCur.subtract(dwBase);
      Map.addLayer(dwDelta.updateMask(dwDelta.abs().gte(DW_CHANGE_MIN)),
        {min: -0.6, max: 0.6, palette: ['#67001f', '#b2182b', '#ef8a62',
                                        '#67a9cf', '#2166ac', '#053061']},
        'DW water probability change (|change| >= ' + DW_CHANGE_MIN + ')',
        false);
    }

    // Agreement between the two definitions is the honest uncertainty band:
    // where they differ is the shallow margin this figure is about.
    if (WATER_METHOD === 'both') {
      var a = pairs[0].dry, b = pairs[1].dry;
      Map.addLayer(a.neq(b).selfMask(), {palette: ['#ff9f0a']},
        YEAR_DRY + ' - NDWI and DW DISAGREE', false);
      print('Area where the two water definitions disagree in ' + YEAR_DRY +
            ' [km2]:', px.updateMask(a.neq(b)).reduceRegion({
        reducer: ee.Reducer.sum(), geometry: aoi, scale: 10,
        maxPixels: 1e9, bestEffort: true, tileScale: 4}).get('area'));
    }

    pairs.forEach(function(w) {
      print(w.tag + ' - water area [km2], baseline then ' + YEAR_DRY + ':',
        ee.Dictionary({
          base_km2: px.updateMask(w.base).reduceRegion({
            reducer: ee.Reducer.sum(), geometry: aoi, scale: 10,
            maxPixels: 1e9, bestEffort: true, tileScale: 4}).get('area'),
          dry_km2: px.updateMask(w.dry).reduceRegion({
            reducer: ee.Reducer.sum(), geometry: aoi, scale: 10,
            maxPixels: 1e9, bestEffort: true, tileScale: 4}).get('area')
        }));
      Map.addLayer(w.base.and(w.dry.not()).selfMask(), {palette: ['#ff3b30']},
        w.tag + ' - water in baseline, dry in ' + YEAR_DRY, false);
    });

    // One extent layer PER METHOD, each named. The previous version drew a
    // single unnamed layer from whichever method happened to be last in the
    // list, so what you were looking at changed with WATER_METHOD without the
    // legend saying so.
    pairs.forEach(function(w) {
      Map.addLayer(w.dry.selfMask(), {palette: ['#2196f3']},
        w.tag + ' - ' + YEAR_DRY + ' water extent', false);
    });

    print('CAVEAT: BOTH definitions are conventions, not measurements. NDWI>0 ' +
          'is a hand-set cut-off; Dynamic World is a classifier whose training ' +
          'data contained few shallow turbid margins like this one. Report the ' +
          'DISAGREEMENT area as the uncertainty. Also note Dynamic World ' +
          'begins 2015-06-27.');
    print('Older note: an NDWI threshold at 0 is a convention, not a measurement. ' +
          'Shallow, turbid and vegetated margins flip either way, and the AOI ' +
          'rectangle clips the lake. Treat the km2 figures as indicative and ' +
          'cite the official Konstanz gauge (BAFU / LUBW) for the actual level.');
  }
}

// ----------------------------------------------------------------------------
// 5b. Legend
// ----------------------------------------------------------------------------
// Only the layers that are actually on the map get an entry, so the legend
// cannot describe something the reader is not looking at.
var legend = ui.Panel({style: {position: 'bottom-left', padding: '8px 12px',
                               backgroundColor: 'white', width: '260px'}});

var head = function(t) {
  return ui.Label({value: t, style: {fontWeight: 'bold', fontSize: '12px',
                                     margin: '6px 0 3px 0'}});
};
// ui.Label has no 'opacity' style property in the Earth Engine UI API, so the
// swatch colour itself carries any transparency.
var row = function(color, text) {
  return ui.Panel({
    widgets: [
      ui.Label({value: '', style: {backgroundColor: color, padding: '7px',
        margin: '0 6px 3px 0', border: '1px solid #999'}}),
      ui.Label({value: text, style: {fontSize: '11px', margin: '0 0 3px 0',
                                     stretch: 'horizontal'}})
    ],
    layout: ui.Panel.Layout.flow('horizontal')
  });
};

legend.add(ui.Label({value: site.name,
  style: {fontWeight: 'bold', fontSize: '13px', margin: '0 0 1px 0'}}));
legend.add(ui.Label({
  value: 'DOY ' + DOY_START + '-' + DOY_END + '  ·  ' + YEAR_DRY +
         ' vs. ' + BASELINE_YEARS[0] + '-' +
         BASELINE_YEARS[BASELINE_YEARS.length - 1] + ' median',
  style: {fontSize: '10px', color: '#555', margin: '0 0 4px 0'}}));

legend.add(head('True colour'));
legend.add(ui.Label({value: 'Identical fixed stretch on both panels — no ' +
  'per-scene levelling, so brightness is comparable between years.',
  style: {fontSize: '10px', color: '#555', margin: '0 0 2px 0'}}));

if (site.water) {
  legend.add(head('Water extent'));
  legend.add(row('#2196f3', YEAR_DRY + ' water'));
  legend.add(row('#ff3b30', 'water in baseline, dry in ' + YEAR_DRY));
  if (WATER_METHOD === 'both') {
    legend.add(row('#ff9f0a', 'the two methods disagree here'));
    legend.add(ui.Label({
      value: 'NDWI>0 is a hand-set cut-off; Dynamic World is a classifier. ' +
             'The orange band is the uncertainty, not an error.',
      style: {fontSize: '10px', color: '#555', margin: '2px 0 2px 0'}}));
  }

  if (WATER_METHOD !== 'ndwi') {
    legend.add(head('DW water probability change'));
    var ramp = ui.Panel({layout: ui.Panel.Layout.flow('horizontal')});
    ['#67001f','#b2182b','#ef8a62','#67a9cf','#2166ac','#053061']
      .forEach(function(c) {
        ramp.add(ui.Label({value: '', style: {backgroundColor: c,
          padding: '7px 16px', margin: '0'}})); });
    legend.add(ramp);
    legend.add(ui.Panel({
      widgets: [
        ui.Label({value: 'less water', style: {fontSize: '10px', margin: '1px 0 0 0'}}),
        ui.Label({value: 'more water', style: {fontSize: '10px',
          margin: '1px 0 0 0', textAlign: 'right', stretch: 'horizontal'}})],
      layout: ui.Panel.Layout.flow('horizontal')}));
    legend.add(ui.Label({
      value: 'Changes below ' + DW_CHANGE_MIN + ' are hidden as model noise.',
      style: {fontSize: '10px', color: '#555', margin: '2px 0 0 0'}}));
  }
}

legend.add(ui.Label({
  value: 'Sentinel-2 L2A · Cloud Score+ cs_cdf >= ' + CLEAR_THRESHOLD +
         ' · baseline includes the 2022 drought, so the anomaly is conservative.',
  style: {fontSize: '9px', color: '#777', margin: '6px 0 0 0'}}));

Map.add(legend);

// ----------------------------------------------------------------------------
// 6. Export — 10 m, LV95
// ----------------------------------------------------------------------------
Export.image.toDrive({
  image: imgDry.visualize(visTrue),
  description: 'S2_' + ACTIVE_SITE + '_' + YEAR_DRY + '_trueColour_LV95',
  folder: 'GEE_Exports', crs: 'EPSG:2056', scale: 10,
  region: aoi, maxPixels: 1e9
});
Export.image.toDrive({
  image: imgBase.visualize(visTrue),
  description: 'S2_' + ACTIVE_SITE + '_baseline_trueColour_LV95',
  folder: 'GEE_Exports', crs: 'EPSG:2056', scale: 10,
  region: aoi, maxPixels: 1e9
});
