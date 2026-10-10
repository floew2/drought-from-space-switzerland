// =============================================================================
// Purity sensitivity — does the answer depend on where we set the threshold?
//
// The question behind it: we keep only cells that are at least PURITY_MIN
// inside a parcel. Why 0.9 and not 0.75, or 0.95? Rather than argue from an
// abstract contamination measure, this measures the effect on the quantity we
// actually report — the NDMI anomaly.
//
// TWO THINGS MAKE THIS CHEAP AND HONEST.
//
// 1. ONE extraction, not one per threshold. The purity band is quantised into
//    classes and a GROUPED reducer returns the mean and pixel count per class
//    for every date in a single pass. The mean at any threshold is then a
//    weighted combination of the classes at or above it, computed in Python.
//    Running the extraction five times would cost five times as much and tell
//    us the same thing.
//
// 2. A FIXED PARCEL SET. Pure cells are not a random subset — they are the
//    interiors of large, compact parcels. So if the anomaly drifts as the
//    threshold falls, that could be contamination, or it could be that small
//    fields simply behave differently. To separate them, only parcels that
//    contain at least one cell above STRICT_PURITY are used, and the threshold
//    is then varied WITHIN that fixed set. Composition stays constant, so any
//    drift is contamination alone.
// =============================================================================

// ----------------------------------------------------------------------------
// 1. Configuration
// ----------------------------------------------------------------------------
var ACTIVE_SITE = 'linthebene';   // 'linthebene' | 'rheintal' | 'alpstein'
var REGIME      = 'intensive';    // matches a mask band, see MASK_BANDS

var YEAR_CUR    = 2026;
var BASELINE_YEARS = [2019, 2020, 2021, 2022, 2023, 2024, 2025];

var DOY_START = 121;   // 1 May
var DOY_END   = 243;   // 31 Aug  — the window the anomaly is quoted over

var INDEX_NAME = 'NDMI';

// Parcels must hold at least one cell this pure to enter the fixed set.
// Under the PSF formulation nothing reaches 1.0, so this cannot be 1.
var STRICT_PURITY = 0.90;

// Purity classes: 0.05-wide bins. Finer than the thresholds anyone would
// actually choose, so the reconstructed curve is not limited by the binning.
var N_CLASSES = 20;

var CLEAR_THRESHOLD = 0.50;
var MAX_SCENE_CLOUD = 80;
var SCALE = 10;   // native mask resolution — no pyramiding of purities

var MASK_ASSET   = 'projects/YOUR_EE_PROJECT/assets/grassland_mask_sg';
var PARCEL_ASSET = 'projects/YOUR_EE_PROJECT/assets/nutzungsflaechen_2025';
var FIELD_CODE   = 'lnf_code';

// Four float bands of PSF-weighted purity 0-1, on the 10 m Sentinel grid.
// Persistence is NOT in the raster — it is in the parcel selection that built
// the mask (grassland in every year AND the same management category), so the
// bands carry only purity and the sensor halo is intact.
var MASK_BANDS = ['mown', 'grazed', 'intensive', 'extensive'];

// Threshold here, not in the mask: sweeping this costs one extraction, not a
// rebuild and re-upload. Measured: the NDMI anomaly varies by 0.5% of its own
// size across the full range 0.05-0.95, at two independent sites.
var PURITY_MIN = 0.9;

var CODES = {
  intensive: [613, 650],
  extensive: [611, 612, 634, 635, 693, 694],
  mown:      [611, 612, 613, 634, 635, 650, 694, 621, 622, 623, 935],
  grazed:    [616, 617, 618, 660, 693, 625, 930, 933]
};

var SITES = {
  linthebene: ee.Geometry.Rectangle([8.90, 47.16, 9.12, 47.30]),
  rheintal:   ee.Geometry.Rectangle([9.50, 47.32, 9.72, 47.50]),
  alpstein:   ee.Geometry.Rectangle([9.25, 47.15, 9.48, 47.30])
};
var aoi = SITES[ACTIVE_SITE];
Map.centerObject(aoi, 12);

// ----------------------------------------------------------------------------
// 2. Fixed parcel set
// ----------------------------------------------------------------------------
var purity = ee.Image(MASK_ASSET).rename(MASK_BANDS)
  .select(REGIME);

var parcels = ee.FeatureCollection(PARCEL_ASSET)
  .filterBounds(aoi)
  .filter(ee.Filter.inList(FIELD_CODE, CODES[REGIME]));

// A parcel qualifies if ANY cell inside it is at least STRICT_PURITY. Those are
// the parcels large and compact enough to have a clean interior — and it is
// exactly within these that the threshold is varied.
var qualified = purity.gte(STRICT_PURITY).reduceRegions({
  collection: parcels,
  reducer: ee.Reducer.max().setOutputs(['pure']),
  scale: SCALE, tileScale: 4
}).filter(ee.Filter.eq('pure', 1));

print('Parcels in AOI:', parcels.size());
print('With at least one cell >= ' + STRICT_PURITY + ':', qualified.size());
print('^ the threshold is varied only WITHIN these, so composition is fixed ' +
      'and any drift is contamination, not a change of sample.');

var eligible = ee.Image(0).byte().paint(qualified, 1).selfMask();
Map.addLayer(qualified, {color: '#3fb950'}, 'fixed parcel set', false);

// Purity class per cell: 0..N_CLASSES-1, restricted to the fixed set.
var pclass = purity.multiply(N_CLASSES).floor()
  .min(N_CLASSES - 1)
  .updateMask(purity.gt(0))
  .updateMask(eligible)
  .rename('pclass').toInt();

// ----------------------------------------------------------------------------
// 3. Imagery
// ----------------------------------------------------------------------------
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED');
var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');

var indexFor = function(img) {
  var m = img.updateMask(img.select('cs_cdf').gte(CLEAR_THRESHOLD));
  return (INDEX_NAME === 'NDVI'
    ? m.normalizedDifference(['B8', 'B4'])
    : m.normalizedDifference(['B8A', 'B11'])).rename('v');
};

var yearCol = function(year) {
  return s2.filterBounds(aoi)
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .filter(ee.Filter.calendarRange(DOY_START, DOY_END, 'day_of_year'))
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', MAX_SCENE_CLOUD))
    .linkCollection(csPlus, ['cs_cdf'])
    .map(indexFor);
};

// ----------------------------------------------------------------------------
// 4. Grouped reduction — one pass gives every threshold
// ----------------------------------------------------------------------------
// Mean index AND pixel count per purity class, per year. The mean at threshold
// t is then sum(mean_c * n_c) / sum(n_c) over classes c >= t, reconstructed in
// Python. That is exactly what a re-extraction at threshold t would return.
var perYear = function(year) {
  var comp = yearCol(year).median().rename('v').updateMask(pclass.mask());
  var stats = comp.addBands(pclass).reduceRegion({
    reducer: ee.Reducer.mean().combine(ee.Reducer.count(), '', true)
      .group({groupField: 1, groupName: 'pclass'}),
    geometry: aoi, scale: SCALE, maxPixels: 1e9, bestEffort: true, tileScale: 4
  });
  return ee.FeatureCollection(
    ee.List(stats.get('groups')).map(function(g) {
      g = ee.Dictionary(g);
      return ee.Feature(null, {
        year: year,
        pclass: g.get('pclass'),
        purity_lo: ee.Number(g.get('pclass')).divide(N_CLASSES),
        mean: g.get('mean'),
        n_px: g.get('count')
      });
    }));
};

var rows = ee.FeatureCollection(
  [YEAR_CUR].concat(BASELINE_YEARS).map(perYear)).flatten();

Export.table.toDrive({
  collection: rows,
  description: 'purity_sensitivity_' + ACTIVE_SITE + '_' + REGIME,
  folder: 'GEE_Exports',
  fileNamePrefix: 'purity_sensitivity_' + ACTIVE_SITE + '_' + REGIME,
  fileFormat: 'CSV',
  selectors: ['year', 'pclass', 'purity_lo', 'mean', 'n_px']
});

print('One task created. The CSV holds mean ' + INDEX_NAME + ' and pixel count '
      + 'per purity class per year; figures.figure_purity_sensitivity() '
      + 'reconstructs the anomaly at every threshold from it.');
