// =============================================================================
// Per-parcel grassland time series — input for CUT DETECTION
//
// Replaces the earlier gee_grassland_extract.js, which predated the mask asset
// and pulled a single index.
//
// WHY PER-PARCEL AND NOT PER-SITE. The trajectory script averages hundreds of
// parcels cut on different dates, which smooths the sawtooth away almost
// completely — that is why mown and grazed curves looked alike there. A cut is
// a parcel-level event and has to be measured at parcel level.
//
// WHAT WE EXPECT TO FIND (2026 Swiss reporting, useful as ground truth):
//   * first and second cuts happened, and early;
//   * the THIRD cut partly failed;
//   * cuts after roughly mid-July did not happen at all;
//   * farms that normally take 5-6 cuts report 2-3.
// So the target signal is the ABSENCE of the usual sawtooth after July, not a
// subtle shift in one event. Missing events are easier to detect than present
// ones, and the normal early cuts act as an internal check that the detector
// works at all.
//
// THIS SCRIPT ONLY EXTRACTS. Detection happens in Python, because the threshold
// will need retuning several times and re-extracting each time is unaffordable.
// =============================================================================

// ----------------------------------------------------------------------------
// 1. Config
// ----------------------------------------------------------------------------
var SITES = {
  linthebene: ee.Geometry.Rectangle([8.90, 47.16, 9.12, 47.30]),
  rheintal:   ee.Geometry.Rectangle([9.50, 47.32, 9.72, 47.50])
};
// Linthebene first: it has by far the best observation density (~51 growing-
// season acquisitions/yr, ~3-day mean gap, against Rheintal's 34 and ~5 days).
// Cut detection needs revisit more than it needs anything else.
var SITE = 'linthebene';

// ORDER MATTERS. Tasks run in submission order, so if quota runs out partway
// you keep whatever finished. Newest first means a partial run still gives
// 2026 plus the most recent baseline years — usable. Oldest first would leave
// you with 2019-2021 and no current year, which is worthless.
var YEARS = [2026, 2025, 2024, 2023, 2022, 2021, 2020, 2019];

// SMOKE TEST FIRST. Per-parcel reduceRegions is heavier than the site-mean
// reduceRegion — check the EECU burn on one year before launching all eight.
var SMOKE_TEST = false;   // full run

// MOWN ONLY. Grazing is continuous defoliation and produces no step to detect;
// including grazed parcels would just add guaranteed false negatives.
var MOWN_CODES = [611, 612, 613, 634, 635, 650, 694,
                  621, 622, 623, 935];

var PARCEL_ASSET = 'projects/YOUR_EE_PROJECT/assets/nutzungsflaechen_2025';
var FIELD_CODE   = 'lnf_code';

// Persistence raster from build_grassland_mask.py. Parcels supply the units,
// the raster restricts which pixels inside them count.
var MASK_ASSET = 'projects/YOUR_EE_PROJECT/assets/grassland_mask_sg';
// Four float bands of PSF-weighted purity 0-1, on the 10 m Sentinel grid.
// Persistence is NOT in the raster — it is in the parcel selection that built
// the mask (grassland in every year AND the same management category), so the
// bands carry only purity and the sensor halo is intact.
var MASK_BANDS = ['mown', 'grazed', 'intensive', 'extensive'];

// Threshold here, not in the mask: sweeping this costs one extraction, not a
// rebuild and re-upload. Measured: the NDMI anomaly varies by 0.5% of its own
// size across the full range 0.05-0.95, at two independent sites.
// mown/grazed is about whether a CUT SIGNAL exists; intensive/extensive is
// about the FARMING SYSTEM, and the two do not nest: 'mown' pools 613 (4-6 cuts
// from early May) with 611/612 (~2 cuts, none before 15 June).

var DOY_START = 60;    // 1 March — need the pre-first-cut buildup
var DOY_END   = 305;   // 1 November

// 250 rather than 400. A mean over 250 parcels is statistically indistinguishable
// from one over 400 for this purpose, and it is ~40% less compute. There is no
// reason to pay for precision the conclusion does not need.
var SAMPLE_N      = 250;
var MIN_AREA_HA   = 0.5;   // before buffering
var EDGE_BUFFER_M = -15;   // mixed edge pixels destroy the step signal
var SCALE = 10;   // native mask resolution — no pyramiding of purities.
                  // A 0.5 ha parcel gives ~50 pixels at this scale.
                          // after buffering, which is plenty for a median, and
                          // this alone cuts pixel reads fourfold.
var SEED          = 42;

// 0.50 rather than 0.60. Catalogue guidance is 0.50-0.65, so this stays within
// range, and it buys roughly 15-25% more dates in the cloudy baseline years —
// which is where the analysis is weakest, since detection sensitivity scales
// with image count. The extra haze that slips through is caught downstream by
// the NDVI-vs-NDMI disagreement test.
var CLEAR_THRESHOLD = 0.50;
var MAX_SCENE_CLOUD = 80;     // drop hopeless scenes before any reduction

// Restrict to ONE MGRS tile. The AOI is covered by two (T32TMT, T32TNT), so
// roughly half the 276 bands are duplicate coverage of the same ground. Set
// this to whichever tile the printed histogram shows carrying more granules;
// parcels outside it drop out of the sample, which is fine for a random sample.
// null = keep both tiles (and pay for both).
var MGRS_TILE = null;   // e.g. 'T32TNT'

// SCALING TEST. The 10-parcel run cost 5,608 EECU-s, but 10 parcels is only
// ~62,000 pixel reads — trivial. So that figure is almost certainly fixed
// overhead (the linkCollection join, building a 276-band stack), which is paid
// whether you ask for 10 parcels or 400. This run tests that:
//   ~6,000  -> overhead-dominated, the full run is affordable
//   ~50,000 -> cost scales with parcels, and the design has to shrink
if (SMOKE_TEST) { YEARS = [2026]; SAMPLE_N = 100; }

// Measured so far, for reference when the tasks land:
//   10 parcels x 1 year, toBands @20m: 5,608 EECU-s (completed, 5 min)
//   100 parcels x 1 year:              >27,900 at 2 min, still climbing
// So cost scales with parcels. This full run is a deliberate gamble on batch
// quota, not a prediction that it is cheap.

var aoi = SITES[SITE];
print('Site:', SITE, '| years:', YEARS, '| sample:', SAMPLE_N);

// ----------------------------------------------------------------------------
// 2. Parcels
// ----------------------------------------------------------------------------
var parcels = ee.FeatureCollection(PARCEL_ASSET)
  .filterBounds(aoi)
  .filter(ee.Filter.inList(FIELD_CODE, MOWN_CODES))
  .map(function(f) {
    return f.set('area_ha', f.geometry().area(1).divide(10000));
  })
  .filter(ee.Filter.gte('area_ha', MIN_AREA_HA))
  .map(function(f) {
    return f.setGeometry(f.geometry().buffer(EDGE_BUFFER_M, 1));
  });

print('Mown parcels >= ' + MIN_AREA_HA + ' ha in AOI:', parcels.size());

// Elevation per parcel: cutting regime is strongly altitude-dependent (roughly
// 5-6 cuts in the lowlands down to 1-2 near the Soemmerung limit), so counts
// cannot be pooled across elevation without stratifying.
var elev = ee.Image('CGIAR/SRTM90_V4').select('elevation');
var withElev = elev.reduceRegions({
  collection: parcels,
  reducer: ee.Reducer.mean().setOutputs(['elev_m']),
  scale: 90, tileScale: 4
}).filter(ee.Filter.notNull(['elev_m']));

var sample = withElev.randomColumn('rnd', SEED).sort('rnd').limit(SAMPLE_N)
  .map(function(f) { return f.set('parcel_id', f.get('system:index')); });

Map.centerObject(aoi, 12);
Map.addLayer(sample, {color: '#ff3b30'}, 'sampled parcels');

// ----------------------------------------------------------------------------
// 3. Imagery
// ----------------------------------------------------------------------------
// The mask band is CONTINUOUS purity 0-1: the fraction of each 20 m cell inside
// a parcel. Threshold here, so a sweep costs one extraction rather than a mask
// rebuild and re-upload.
// NOTE: GEE's default pyramiding for float is 'mean', so reading this at a
// coarser scale would average purities — a plausible number that is the purity
// of nothing. Upload with pyramidingPolicy {'.default': 'sample'}, or reduce at
// the native 20 m.
var PURITY_MIN = 0.9;
var persistMask = ee.Image(MASK_ASSET).rename(MASK_BANDS)
  .select('mown').gte(PURITY_MIN);

var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED');
var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');

// NDVI carries the cut signal; NDMI is the HAZE DISCRIMINATOR and is the reason
// both are exported. Residual thin haze depresses NDVI while barely touching
// NDMI — at Rheintal in 2026, NDVI fell 0.197 in two days while NDMI moved
// 0.028. That looks exactly like a cut to any NDVI-threshold detector, so it
// would generate false positives. A REAL cut removes biomass and drops BOTH.
// The disagreement between the two is therefore the test, and it is something
// most published cut-detection work does not do.
// COST NOTE — two failed attempts, and what they taught.
//   1st: one reduceRegions per granule. 60 parcels x 1 year = ~8,500 EECU-s.
//   2nd: mosaicking granules by date first. WORSE (~31,800), because
//        ee.Join + mosaic() builds a COMPUTED image that GEE does not cache,
//        so every reduction re-derived it from its granules.
// The parcels total ~1 km2, so ~3M pixel reads per year — that is nearly free.
// The cost was therefore per-CALL overhead on ~17,700 separate reduceRegions
// invocations, not pixel work. So: stack the whole season into ONE image with
// toBands() and make ONE reduceRegions call per year.
//
// Also dropped kNDVI and NDRE from the export. kNDVI = tanh(NDVI^2) is
// derivable in Python at zero cost, and the detector never uses NDRE — so
// exporting them was doubling the band count for nothing.
var indexCol = function(year) {
  var raw = s2
    .filterBounds(aoi)
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .filter(ee.Filter.calendarRange(DOY_START, DOY_END, 'day_of_year'))
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', MAX_SCENE_CLOUD));
  if (MGRS_TILE) {
    raw = raw.filter(ee.Filter.eq('MGRS_TILE', MGRS_TILE.replace('T', '')));
  }
  raw = raw.linkCollection(csPlus, ['cs_cdf']);

  return raw.map(function(img) {
      var m = img.updateMask(img.select('cs_cdf').gte(CLEAR_THRESHOLD))
                 .updateMask(persistMask);
      var ndvi = m.normalizedDifference(['B8',  'B4']).rename('NDVI');
      var ndmi = m.normalizedDifference(['B8A', 'B11']).rename('NDMI');
      // Band names must survive toBands(), so the date and tile go into
      // system:index. Result: "20250315_32TNT_NDVI".
      var tag = ee.Image(img).date().format('YYYYMMdd')
        .cat('_').cat(ee.String(ee.Image(img).get('MGRS_TILE')));
      return ndvi.addBands(ndmi).set('system:index', tag);
    });
};

YEARS.forEach(function(y) {
  print('Granules, ' + y + ':', indexCol(y).size());
});
// Which MGRS tiles cover this AOI? If one tile covers it adequately, setting
// MGRS_TILE removes duplicate coverage at zero cost.
print('MGRS tiles present:', s2.filterBounds(aoi)
  .filter(ee.Filter.calendarRange(YEARS[0], YEARS[0], 'year'))
  .aggregate_histogram('MGRS_TILE'));

// ----------------------------------------------------------------------------
// 4. Export — one task per year
// ----------------------------------------------------------------------------
// Median over each parcel's interior, not mean: robust to a stray bright pixel.
// A band is empty for a parcel-date when cloud left no valid pixels there, and
// the Python side simply drops those — a PER-PARCEL judgement, since a 2 ha
// field does not care whether the wider site was clear.
YEARS.forEach(function(year) {
  // ONE image, ONE reduceRegions. Bands are <date>_<tile>_<index>.
  var stack = indexCol(year).toBands();

  // median AND count. The count was missing from earlier exports, which left
  // MIN_PARCEL_COVERAGE in the Python code doing nothing at all: with no pixel
  // count, a parcel-date was kept whenever ANY pixel was valid, so heavily
  // clouded parcels contributed on equal terms with clear ones. Doubling the
  // columns is cheap; both statistics come from the same pass.
  var rows = stack.reduceRegions({
    collection: sample,
    reducer: ee.Reducer.median().combine(ee.Reducer.count(), '', true),
    scale: SCALE,
    tileScale: 8
  });

  Export.table.toDrive({
    collection: rows,
    description: 'cuts_' + SITE + '_' + year + (SMOKE_TEST ? '_SMOKE' : ''),
    folder: 'GEE_Exports',
    fileNamePrefix: 'cuts_' + SITE + '_' + year + (SMOKE_TEST ? '_SMOKE' : ''),
    fileFormat: 'CSV'
    // No selectors: band names are generated, so let everything through and
    // reshape wide -> long in Python.
  });
});

print('One task per year. Output is WIDE: one row per parcel, and now TWO ' +
      'columns per date/tile/index — <name>_median and <name>_count. ' +
      'reshape_wide() handles both this and the older single-column format.');
print('RUN THE SMOKE TEST FIRST and watch the EECU counter for two minutes. ' +
      'Two cost predictions have already been wrong; measure, do not trust.');
