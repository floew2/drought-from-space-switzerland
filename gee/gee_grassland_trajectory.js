// =============================================================================
// Grassland VI trajectory extraction — site-level seasonal curves
//
// Reads the mask asset built by drought.mask.build_mask() and uploaded by hand. Four float
// bands of PSF-weighted purity 0-1 on the 10 m Sentinel grid; persistence is
// already in the parcel selection, so the raster carries purity alone.
//
// Produces ONE small CSV: mean index over the grassland mask, per acquisition
// date, per year, per site, split by regime (mown / grazed).
//
// COST NOTE — an early version burned ~97,000 EECU-seconds in four minutes and
// had to be killed, because the polygon-to-raster conversion sat INSIDE the
// per-image loop and was redone for every image/regime/year combination. The
// mask is now a stored asset and is simply read.
//
// Downstream: plot in Python with the same styling as the precipitation figure
// (trajectory panel + cumulative-integral panel).
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
var CONFIG_HASH = "3253fe2577fa";   // settings + mask this copy was rendered with; goes into export names
var PURITY_MIN = 0.9;   // mask.purity_min
var CLEAR_THRESHOLD = 0.6;   // Cloud Score+ cs_cdf, gee.clear_threshold
var MASK_BANDS = ["mown", "grazed", "intensive", "extensive"];   // band order of the mask raster, as written by build_mask()
var REGIMES = ["intensive", "extensive"];   // analysis.regimes
var YEARS = [2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026];   // project.baseline_years + current_year
// <<< end of generated block

var SITES = {
  rheintal: {
    name: 'St. Galler Rheintal',
    geom: ee.Geometry.Rectangle([9.50, 47.32, 9.72, 47.50]),
    stratum: 'lowland'
  },
  linthebene: {
    name: 'Linthebene',
    geom: ee.Geometry.Rectangle([8.90, 47.16, 9.12, 47.30]),
    stratum: 'lowland'
  },
  alpstein: {
    name: 'Alpstein / Churfirsten',
    geom: ee.Geometry.Rectangle([9.25, 47.15, 9.48, 47.30]),
    stratum: 'summering'
  }
};

// Run all three in one go — it is cheap enough, and one CSV keeps the Python
// side simple. Drop entries here if you want to shorten it.
var ACTIVE_SITES = ['rheintal', 'linthebene', 'alpstein'];

// SMOKE TEST FIRST. Set true, run, and check the printed row count and the
// EECU burn in the Tasks tab before launching the full job. One site, two
// years. If this is cheap, the full run is roughly 12x it.
var SMOKE_TEST = false;   // smoke test passed at 338 EECU-s; full run ~4,000

// YEARS: rendered from config.yaml (baseline + current year).

if (SMOKE_TEST) {
  ACTIVE_SITES = ['rheintal'];
  YEARS = [2025, 2026];
}

// ONE mask asset, built locally by drought.mask.build_mask() and uploaded.
// GeoTIFF band descriptions are not reliably carried into GEE, so the bands are
// renamed explicitly here BY POSITION. Order must match the writer:
//   b1 mown | b2 grazed | b3 intensive | b4 extensive
// MASK_ASSET: rendered from config.yaml.
// Band order must match what drought.mask.build_mask() printed — GEE names them b1..bN, so
// they are renamed BY POSITION. Do not reorder.
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

// CONTINUOUS purity 0-1: the share of each 10 m cell's sensor RESPONSE that
// comes from inside a parcel, PSF-weighted (sigma 5.3 m from the specified
// system MTF). Threshold here rather than in the mask, so a sweep costs one
// cheap extraction each instead of a rebuild and re-upload.
// PURITY_MIN: rendered from config.yaml (mask.purity_min).
// Upload the asset with pyramidingPolicy {'.default': 'sample'}: the default
// for float is 'mean', which would average purities at coarser scales.

// Persistence-filtered mask: grassland in EVERY year, not just the regime year.
// Stops a parcel that was arable in 2019 contributing an arable curve to the
// grassland baseline. Set false to quantify how much the choice matters.

// FULL season, not the anomaly window: the trajectory panel exists to show WHEN
// 2026 diverged, which needs the curve either side of the divergence. The
// alpine site also needs the spring shoulder to reveal snowmelt timing.
var DOY_START = 60;    // 1 March
var DOY_END   = 305;   // 1 November

var REDUCE_SCALE  = 10;   // NATIVE mask resolution. Reducing coarser makes GEE
                          // pyramid the float purity band, and an averaged
                          // purity is the purity of nothing. The trajectory ran
                          // at ~4,000 EECU-s, so the finer read is affordable.
var MIN_PIXELS    = 2000; // reject dates where cloud left too little of the
                          // site. Raised with the finer scale: at 10 m a site
                          // holds ~25x the pixels it did at 50 m.

// ----------------------------------------------------------------------------
// 2. Grassland mask — READ, not computed
// ----------------------------------------------------------------------------
// The whole optimisation lives here. ee.Image() on a stored asset reads
// finished pixels; nothing is re-derived from polygons per image. The earlier
// version rasterised parcels inside the per-image loop and burned ~137,000
// EECU-seconds; this one ran the same work for ~170 per site-year.
var maskImg = ee.Image(MASK_ASSET).rename(MASK_BANDS);

// Which grouping to extract. 'mown'/'grazed' is about whether a cut signal
// exists; 'intensive'/'extensive' is about the farming system. Use intensive
// for the cut work — pooling 613 with 611/612 averages two systems with very
// different cutting schedules.
// REGIMES: rendered from config.yaml (analysis.regimes).

var maskFor = function(regime) {
  return maskImg
    .select(regime)
    .gte(PURITY_MIN)
    .selfMask();
};

print('Mask asset bands:', maskImg.bandNames());
print('| purity threshold:', PURITY_MIN);

// ----------------------------------------------------------------------------
// 3. Sentinel-2 index collection
// ----------------------------------------------------------------------------
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED');
var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');

var indexCol = function(geom, year) {
  return s2
    .filterBounds(geom)
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .filter(ee.Filter.calendarRange(DOY_START, DOY_END, 'day_of_year'))
    .linkCollection(csPlus, ['cs_cdf'])
    .map(function(img) {
      var m = img.updateMask(img.select('cs_cdf').gte(CLEAR_THRESHOLD));
      var ndvi = m.normalizedDifference(['B8',  'B4']).rename('NDVI');
      var ndmi = m.normalizedDifference(['B8A', 'B11']).rename('NDMI');
      var ndre = m.normalizedDifference(['B8A', 'B5']).rename('NDRE');
      // kNDVI with sigma = 0.5*(NIR+Red) reduces to tanh(NDVI^2).
      var kndvi = ndvi.pow(2).tanh().rename('kNDVI');
      return ndvi.addBands([kndvi, ndmi, ndre])
        .copyProperties(img, ['system:time_start']);
    });
};

// ----------------------------------------------------------------------------
// 4. Build the table: one row per site x regime x year x date
// ----------------------------------------------------------------------------
var rows = ee.FeatureCollection(ACTIVE_SITES.map(function(siteKey) {
  var site = SITES[siteKey];

  var perRegime = REGIMES.map(function(regime) {
    var mask = maskFor(regime);

    var perYear = YEARS.map(function(year) {
      return indexCol(site.geom, year).map(function(img) {
        var masked = ee.Image(img).updateMask(mask);
        // ONE reduceRegion returns the mean of all four bands plus a count.
        // Four indices therefore cost essentially the same as one.
        var st = masked.reduceRegion({
          reducer: ee.Reducer.mean().combine(ee.Reducer.count(), '', true),
          geometry: site.geom,
          scale: REDUCE_SCALE,
          maxPixels: 1e8,
          bestEffort: true,
          tileScale: 4
        });
        return ee.Feature(null, {
          site: siteKey,
          regime: regime,
          year: year,
          date: ee.Image(img).date().format('YYYY-MM-dd'),
          doy: ee.Image(img).date().getRelative('day', 'year').add(1),
          NDVI:  st.get('NDVI_mean'),
          kNDVI: st.get('kNDVI_mean'),
          NDMI:  st.get('NDMI_mean'),
          NDRE:  st.get('NDRE_mean'),
          // count from NDVI: all four share the same cloud mask, so one count
          // describes them all.
          n_px:  st.get('NDVI_count')
        });
      });
    });
    return ee.FeatureCollection(perYear).flatten();
  });
  return ee.FeatureCollection(perRegime).flatten();
})).flatten()
  .filter(ee.Filter.notNull(['NDVI', 'kNDVI', 'NDMI', 'NDRE']))
  .filter(ee.Filter.gte('n_px', MIN_PIXELS));

// ----------------------------------------------------------------------------
// 5. Export
// ----------------------------------------------------------------------------
Export.table.toDrive({
  collection: rows,
  description: 'grassland_trajectory_multi_p' +
               Math.round(PURITY_MIN * 100) + (SMOKE_TEST ? '_SMOKE' : ''),
  folder: 'GEE_Exports',
  fileNamePrefix: 'grassland_trajectory_multi_p' +
                  Math.round(PURITY_MIN * 100) + (SMOKE_TEST ? '_SMOKE' : '') +
                  '_cfg' + CONFIG_HASH,
  fileFormat: 'CSV',
  selectors: ['site', 'regime', 'year', 'date', 'doy',
              'NDVI', 'kNDVI', 'NDMI', 'NDRE', 'n_px']
});

print('SMOKE_TEST =', SMOKE_TEST, '| sites:', ACTIVE_SITES, '| years:', YEARS);
print('One export task created. Expect a few thousand rows — small enough to ' +
      'open in a text editor and sanity-check before plotting.');
print('Clear-observation counts per year come out of the CSV itself (count ' +
      'rows per site/regime/year). Report that table alongside the curves: if ' +
      '2026 has MORE clear dates than the baseline years — which sunny drought ' +
      'summers produce — then any deficit was measured under better observing ' +
      'conditions, not worse.');
