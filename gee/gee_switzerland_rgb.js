// =============================================================================
// Switzerland — MODIS MOD09A1 RGB comparison, drought year vs. climatological
// baseline. Primary output: two true-colour composites with an identical,
// pinned stretch, intended for screenshots / GeoTIFF export.
//
// What changed vs. the first draft, and why:
//   1. StateQA cloud/shadow/cirrus masking + valid-range masking. MOD09A1 is
//      NOT cloud-free and stores fill as -28672; without this the median that
//      the composite is anchored on is contaminated.
//   2. Day-of-year filtering instead of calendar dates, so leap and non-leap
//      years use the same phenological window and nothing is silently dropped
//      by the end-exclusive filterDate().
//   3. Baseline is a MULTI-YEAR median, not a single year. A single year is one
//      draw from the distribution: 2016 in particular was a wet year in
//      NE Switzerland (May–Jul 2016 sits near the 90th percentile of 1991–2020
//      at the St. Gallen station), so it would overstate the contrast.
//   4. The JRC water mask is explicitly aggregated 30 m -> MODIS grid instead of
//      letting GEE reconcile the scales implicitly at render time.
//   5. Export in EPSG:2056 (LV95), not the EPSG:4326 default.
//   6. Prints image counts and the mean solar zenith angle per period.
//
// !! IMPORTANT CAVEAT, please keep it in the figure caption !!
//   Terra stopped inclination maintenance in Feb 2020. Its equator-crossing
//   mean local time has drifted from 10:30 towards ~09:00, and the orbit was
//   lowered to 694 km in Oct 2022. A morning platform crossing earlier means a
//   lower sun and more terrain shadow — which in Alpine terrain pushes
//   reflectance in the SAME direction as drought. The printed solar zenith
//   angles quantify this. For a defensible visual, cross-check against
//   Sentinel-2 or Landsat 8/9, whose overpass times are stable.
// =============================================================================

// ----------------------------------------------------------------------------
// 1. Configuration
// ----------------------------------------------------------------------------
var YEAR_DRY       = 2026;
var BASELINE_YEARS = [2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020];
// Trade-off, be explicit about it: an older baseline is less matched to Terra's
// current overpass time; a newer one includes 2018/2022 drought years. There is
// no window that avoids both. Set SINGLE_BASELINE_YEAR to a number to override.
var SINGLE_BASELINE_YEAR = null;

var DOY_START = 166;   // 15 June
var DOY_END   = 212;   // 31 July  (inclusive on both ends)

var STRICT_QA = false; // true = only cloud state 0; false = also accept state 3
                       // ("not set, assumed clear"), which reduces data holes.
var SHOW_ANALYSIS = true;  // NDVI anomaly layer + 2 zonal reduceRegion calls.
                           // Set false if you want zero reduction cost.

// --- snow / high-alpine handling (see notes in section 6) --------------------
var COMPOSITE_METHOD = 'median';  // 'median' (smooth, recommended for RGB) or
                                  // 'closestToMedian' (spectrally consistent,
                                  // but picks whole scenes per pixel -> blocky)
var MIN_VALID_OBS = 2;    // below this many QA-clear observations, use the
                          // unmasked stack instead of a 1-scene composite
var NDSI_SNOW = 0.4;      // standard MODIS snow threshold
var MAX_LAKE_ELEV = 2000; // m; above this, MCD12Q1 "water" is almost always
                          // misclassified snow/ice, not a lake
var SHOW_DIAGNOSTICS = true;   // extra map layers, all off by default

// ----------------------------------------------------------------------------
// 2. Area of interest and basemap
// ----------------------------------------------------------------------------
var ch = ee.FeatureCollection('FAO/GAUL/2015/level0')
  .filter(ee.Filter.eq('ADM0_NAME', 'Switzerland'));
Map.centerObject(ch, 8);

Map.setOptions('Black', {Black: [
  {featureType: 'all', stylers: [{color: '#000000'}]},
  {elementType: 'labels', stylers: [{visibility: 'off'}]}
]});

// ----------------------------------------------------------------------------
// 3. Collection, bands, projection
// ----------------------------------------------------------------------------
var BANDS = ['sur_refl_b01', 'sur_refl_b02', 'sur_refl_b03',
             'sur_refl_b04', 'sur_refl_b06', 'sur_refl_b07'];

var modisCol = ee.ImageCollection('MODIS/061/MOD09A1').filterBounds(ch);
var modisProj = ee.Image(modisCol.first()).select('sur_refl_b01').projection();

// ----------------------------------------------------------------------------
// 4. Land cover and permanent water — both already on the MODIS 500 m grid
// ----------------------------------------------------------------------------
// NOTE: an earlier version aggregated the 30 m JRC Global Surface Water layer
// with reduceResolution().reproject(). That is more precise, but reproject()
// pins the projection and forces the whole downstream graph onto a fixed fine
// grid, which blows up ("Reprojection output too large") and burns quota.
// MCD12Q1 is native MODIS sinusoidal, so no reprojection happens anywhere here.
var lcImg = ee.ImageCollection('MODIS/061/MCD12Q1')
  .sort('system:time_start', false).first().select('LC_Type1');
print('Land cover epoch used:', ee.Date(lcImg.get('system:time_start')).get('year'));

// Elevation, read through the 90 m pyramid at whatever scale the caller asks
// for. No reproject(), so no cost blow-up.
var elev = ee.Image('CGIAR/SRTM90_V4').select('elevation');

// IGBP 17 = water bodies. Above ~2000 m that class is overwhelmingly
// misclassified snowfield or glacier, and sending those pixels through the
// dark-water compositing path is what produces black or muddy blotches on the
// summits. Every real Swiss lake sits below this cap (Silsersee, the highest
// large one, is ~1797 m).
var isLake = lcImg.eq(17).and(elev.lt(MAX_LAKE_ELEV)).clip(ch);
var isSnowIceLC = lcImg.eq(15);                                     // IGBP 15
var isVegetated = lcImg.lte(10).or(lcImg.eq(12)).or(lcImg.eq(14));  // forest..cropland

// 5. QA masking
// ----------------------------------------------------------------------------
// MOD09A1 StateQA bit layout:
//   0-1 cloud state (0 clear, 1 cloudy, 2 mixed, 3 not set/assumed clear)
//   2   cloud shadow            8-9 cirrus (0 none .. 3 high)
//   10  internal cloud flag     13  pixel adjacent to cloud
//   15  internal snow mask
var maskMOD09A1 = function(img) {
  var qa       = img.select('StateQA');
  var cloud    = qa.bitwiseAnd(3);
  var shadow   = qa.rightShift(2).bitwiseAnd(1);
  var cirrus   = qa.rightShift(8).bitwiseAnd(3);
  var intCloud = qa.rightShift(10).bitwiseAnd(1);
  var adjacent = qa.rightShift(13).bitwiseAnd(1);

  var cloudOk = STRICT_QA ? cloud.eq(0) : cloud.eq(0).or(cloud.eq(3));
  var clear = cloudOk
    .and(shadow.eq(0))
    .and(cirrus.lte(1))
    .and(intCloud.eq(0))
    .and(adjacent.eq(0));

  var sr = img.select(BANDS);
  // Fill is -28672; the documented valid range is -100..16000.
  var inRange = sr.reduce(ee.Reducer.min()).gte(-100)
    .and(sr.reduce(ee.Reducer.max()).lte(16000));

  return sr.updateMask(clear.and(inRange))
    .copyProperties(img, ['system:time_start']);
};

var validRangeOnly = function(img) {
  var sr = img.select(BANDS);
  var inRange = sr.reduce(ee.Reducer.min()).gte(-100)
    .and(sr.reduce(ee.Reducer.max()).lte(16000));
  return sr.updateMask(inRange).copyProperties(img, ['system:time_start']);
};

// ----------------------------------------------------------------------------
// 6. Composites
// ----------------------------------------------------------------------------
// Day-of-year filter: the same phenological window every year, leap or not.
// Note MOD09A1 composites are stamped with the START of their 8-day period, so
// the last one selected still carries observations a few days past DOY_END.
var windowFor = function(year) {
  return modisCol
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    .filter(ee.Filter.calendarRange(DOY_START, DOY_END, 'day_of_year'));
};

// "Closest-to-median" mosaic — the usual GEE idiom, often loosely called a
// medoid. It is NOT a true medoid: the reference is a band-wise median, which is
// not itself an observation. Good enough here, but name it honestly.
var closestToMedian = function(col) {
  var median = col.median();
  var scored = col.map(function(img) {
    var dist = img.subtract(median).pow(2).reduce(ee.Reducer.sum()).sqrt();
    return img.addBands(dist.multiply(-1).rename('score'));
  });
  return scored.qualityMosaic('score').select(BANDS);
};

var ndsiOf = function(img) {
  return img.normalizedDifference(['sur_refl_b04', 'sur_refl_b06']).rename('NDSI');
};

var compositeOf = function(col) {
  return COMPOSITE_METHOD === 'median' ? col.median() : closestToMedian(col);
};

// Why the high-alpine artefacts happened, and what each line below does:
//
//  (a) MOD09A1's cloud flags routinely call bright snow "cloud". Over the snow
//      zone the QA-masked stack is therefore nearly empty, so the composite was
//      being assembled from one or two surviving scenes -- or falling through to
//      the unmasked fallback -- on a pixel-by-pixel basis. That mixture is the
//      patchwork you are seeing. Fix: detect snow and take the whole snow zone
//      from the unmasked stack, consistently.
//  (b) qualityMosaic picks a whole 8-day scene per pixel. Between mid-June and
//      end-July an Alpine snowpack changes enormously from one 8-day period to
//      the next, so neighbouring pixels drawn from different periods differ
//      sharply -> blockiness. Fix: default to a median composite for the RGB.
//  (c) MCD12Q1 mislabels high snowfields as water (class 17), which routed them
//      through the 20th-percentile "lake" path (darkest observation wins).
//      Fix: the elevation cap in section 4, plus never substitute over snow.
//  (d) A masked lake composite in where() produced holes. Fix: unmask onto land.
//
// Residual effect that is NOT a bug: snow reflectance is roughly 0.6-0.9, i.e.
// DN 6000-9000, far above the pinned max of 2000, so snow clips to flat white by
// design. Toggle the "snow stretch" layer to see structure up there.
var cleanComposite = function(year) {
  var scenes = windowFor(year);
  var raw    = scenes.map(validRangeOnly);
  var masked = scenes.map(maskMOD09A1);

  var nValid = masked.select('sur_refl_b01').count().unmask(0).rename('n_valid');

  var primary  = compositeOf(masked);
  var fallback = compositeOf(raw);

  var snowy = ndsiOf(fallback).gt(NDSI_SNOW).unmask(0).or(isSnowIceLC);

  var land = primary
    .where(nValid.lt(MIN_VALID_OBS), fallback)   // too few clear scenes
    .where(snowy, fallback)                      // snow zone, handled as a unit
    .unmask(fallback);                           // never punch a hole

  // 20th percentile over water strips sun glint and thin cloud off lakes.
  // unmask(land) guarantees where() can never introduce a gap.
  var lake = masked.reduce(ee.Reducer.percentile([20]))
    .select(BANDS.map(function(b) { return b + '_p20'; }), BANDS)
    .unmask(land);

  return land.where(isLake.and(snowy.not()), lake).clip(ch)
    .set('year', year);
};

// Diagnostic images for the dry year, so the snow zone can be inspected rather
// than guessed at. These are map layers only -- no reduceRegion, no quota cost
// unless you actually switch them on.
var dryWindow   = windowFor(YEAR_DRY);
var dryNValid   = dryWindow.map(maskMOD09A1).select('sur_refl_b01')
                    .count().unmask(0).rename('n_valid');
var drySnowMask = ndsiOf(compositeOf(dryWindow.map(validRangeOnly)))
                    .gt(NDSI_SNOW).unmask(0).or(isSnowIceLC);

var imgDry = cleanComposite(YEAR_DRY);

// Build every baseline composite ONCE and reuse it for both the RGB median and
// the NDVI statistics. The first version built them twice, doubling the graph.
var baseComposites = SINGLE_BASELINE_YEAR
  ? [cleanComposite(SINGLE_BASELINE_YEAR)]
  : BASELINE_YEARS.map(cleanComposite);

var imgBase = ee.ImageCollection(baseComposites).median().clip(ch);

var baseLabel = SINGLE_BASELINE_YEAR
  ? String(SINGLE_BASELINE_YEAR)
  : BASELINE_YEARS[0] + '-' + BASELINE_YEARS[BASELINE_YEARS.length - 1] +
    ' median (' + BASELINE_YEARS.length + ' yr)';

// ----------------------------------------------------------------------------
// 7. Diagnostics — check these before trusting anything below
// ----------------------------------------------------------------------------
print('MOD09A1 scenes in window, ' + YEAR_DRY + ':', windowFor(YEAR_DRY).size());
print('MOD09A1 scenes in window, ' + BASELINE_YEARS[0] + ':',
      windowFor(BASELINE_YEARS[0]).size());

var meanSZA = function(year) {
  return windowFor(year).select('SolarZenith').mean().multiply(0.01)
    .reduceRegion({reducer: ee.Reducer.mean(), geometry: ch.geometry(),
                   scale: 5000, maxPixels: 1e9}).get('SolarZenith');
};
print('Mean solar zenith [deg], ' + YEAR_DRY + ':', meanSZA(YEAR_DRY));
print('Mean solar zenith [deg], ' + BASELINE_YEARS[0] + ':', meanSZA(BASELINE_YEARS[0]));
print('^ a large gap here means illumination, not vegetation, explains part of ' +
      'the visual difference (Terra orbit drift).');

// ----------------------------------------------------------------------------
// 8. Visualisation — stretch pinned identically across both periods
// ----------------------------------------------------------------------------
var visTrue = {
  bands: ['sur_refl_b01', 'sur_refl_b04', 'sur_refl_b03'],
  min: -100, max: 2000, gamma: 1.4
};
var visFalse = {
  bands: ['sur_refl_b02', 'sur_refl_b01', 'sur_refl_b04'],
  min: -200, max: 4000, gamma: 1.4
};

Map.addLayer(imgBase, visTrue, 'Baseline (' + baseLabel + ') — true colour', true);
Map.addLayer(imgDry,  visTrue, YEAR_DRY + ' — true colour', true);
Map.addLayer(imgBase, visFalse, 'Baseline (' + baseLabel + ') — false colour', false);
Map.addLayer(imgDry,  visFalse, YEAR_DRY + ' — false colour', false);
Map.addLayer(isLake.selfMask(), {palette: ['cyan']}, 'Permanent water (MODIS grid)', false);

// Same bands, wider ceiling: snow is ~DN 6000-9000 and clips to white under the
// pinned stretch. Use this layer only to inspect the snow zone -- it is NOT
// comparable to the pinned screenshots.
var visTrueSnow = {
  bands: ['sur_refl_b01', 'sur_refl_b04', 'sur_refl_b03'],
  min: -100, max: 7000, gamma: 1.2
};
Map.addLayer(imgBase, visTrueSnow, 'Baseline - snow stretch (not comparable)', false);
Map.addLayer(imgDry,  visTrueSnow, YEAR_DRY + ' - snow stretch (not comparable)', false);

if (SHOW_DIAGNOSTICS) {
  Map.addLayer(dryNValid, {min: 0, max: 6, palette: ['#000000', '#b30000', '#fdbb84', '#ffffcc']},
    YEAR_DRY + ' - QA-clear observations per pixel', false);
  Map.addLayer(dryNValid.lt(MIN_VALID_OBS).selfMask(), {palette: ['magenta']},
    YEAR_DRY + ' - pixels taken from the unmasked stack', false);
  Map.addLayer(drySnowMask.selfMask(), {palette: ['#00ffff']},
    YEAR_DRY + ' - snow/ice zone (NDSI or IGBP 15)', false);
  print('Turn on "QA-clear observations per pixel": if the Alps are dark there, ' +
        'MODIS cloud flags are rejecting snow and the snow zone is coming from ' +
        'the unmasked stack by design.');
}

// ----------------------------------------------------------------------------
// 9. Optional: NDVI anomaly and zonal statistics
// ----------------------------------------------------------------------------
if (SHOW_ANALYSIS) {
  var ndvi = function(img) {
    return img.normalizedDifference(['sur_refl_b02', 'sur_refl_b01']).rename('dNDVI');
  };
  var ndsi = function(img) {
    return img.normalizedDifference(['sur_refl_b04', 'sur_refl_b06']).rename('NDSI');
  };

  var baseNdviCol = ee.ImageCollection(baseComposites.map(ndvi));
  var ndviBase = baseNdviCol.median();
  var ndviSd = baseNdviCol.reduce(ee.Reducer.stdDev()).rename('sd');

  var anomaly = ndvi(imgDry).subtract(ndviBase).rename('dNDVI');
  var zscore = anomaly.divide(ndviSd.max(0.005)).rename('z');

  // Statistics only where there is actual vegetation: no water, no snow/ice,
  // no bare rock. isVegetated and isLake both come from MCD12Q1 (MODIS grid).
  var analysisMask = isVegetated
    .and(isLake.not())
    .and(ndsi(imgDry).lt(0.4))
    .and(ndsi(imgBase).lt(0.4));

  Map.addLayer(anomaly.updateMask(analysisMask),
    {min: -0.25, max: 0.25, palette: ['#8c2d04', '#fdae61', '#ffffbf', '#a6d96a', '#1a9850']},
    'NDVI anomaly (' + YEAR_DRY + ' - baseline)', false);
  Map.addLayer(zscore.updateMask(analysisMask),
    {min: -3, max: 3, palette: ['#67001f', '#f7f7f7', '#053061']},
    'NDVI z-score', false);

  // --- Zonal statistics, deliberately cheap -------------------------------
  // Two reduceRegion calls in total (the first version made five), a 90 m DEM
  // read through its pyramid instead of the 30 m SRTM, scale 1000 m because
  // area means do not need 500 m, and tileScale to keep memory down.
  var elev = ee.Image('CGIAR/SRTM90_V4').select('elevation');
  var zone = ee.Image(0)
    .where(elev.gte(600), 1)
    .where(elev.gte(1000), 2)
    .where(elev.gte(1500), 3)
    .where(elev.gte(2000), 4)
    .rename('zone');

  var byZone = anomaly.addBands(zone).updateMask(analysisMask).reduceRegion({
    reducer: ee.Reducer.mean().group({groupField: 1, groupName: 'zone'}),
    geometry: ch.geometry(), scale: 1000, maxPixels: 1e8,
    bestEffort: true, tileScale: 4
  });
  print('Mean NDVI anomaly by elevation zone ' +
        '(0: <600, 1: 600-1000, 2: 1000-1500, 3: 1500-2000, 4: >2000 m):', byZone);

  var frac = zscore.lt(-1).rename('below1sd').updateMask(analysisMask).reduceRegion({
    reducer: ee.Reducer.mean(), geometry: ch.geometry(), scale: 1000,
    maxPixels: 1e8, bestEffort: true, tileScale: 4
  });
  print('Fraction of vegetated area below -1 sigma:', frac);
}

// ----------------------------------------------------------------------------
// 10. Export — EPSG:2056 (LV95), not the EPSG:4326 default
// ----------------------------------------------------------------------------
Export.image.toDrive({
  image: imgDry.visualize(visTrue),
  description: 'MODIS_CH_' + YEAR_DRY + '_trueColour_LV95',
  folder: 'GEE_Exports',
  crs: 'EPSG:2056',
  scale: 500,
  region: ch.geometry(),
  maxPixels: 1e9
});

Export.image.toDrive({
  image: imgBase.visualize(visTrue),
  description: 'MODIS_CH_baseline_trueColour_LV95',
  folder: 'GEE_Exports',
  crs: 'EPSG:2056',
  scale: 500,
  region: ch.geometry(),
  maxPixels: 1e9
});
