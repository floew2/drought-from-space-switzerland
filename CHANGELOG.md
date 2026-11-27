# Changelog

## 2.0.1

### Changed
- README, METHOD and gallery report the swissEO run of 27 November 2026, with
  its figures; Alpstein now has season integrals for both meadow types.

## 2.0.0

### Changed
- Satellite data now come from swisstopo's swissEO S2-SR v200, downloaded
  once and computed locally. Google Earth Engine is no longer used; the
  results quoted in the README and METHOD are from the first run in Earth
  Engine and will be replaced by the swissEO run.
- The package is split into `core`, `sources`, `mask`, `products`,
  `analysis`, `viz` and `pipeline`.
- One command runs the pipeline: `drought run [STAGE]` runs a stage and
  whatever it needs that is not done, judged by each output's provenance record.
- Commands: `drought fetch`, `drought mask build`, `drought condition`,
  `sensitivity`, `anomaly`, `cuts`, `vhi`, `evaluate`, `precip`, `report`,
  `inspect`. The `swisseo`, `ee`, `trajectory` and `parcels` commands are gone.
- Data folders: `data/scenes/<site>/` (downloaded mosaics) and
  `data/series/<site>/` (tables). `trajectory.csv` is now `condition.csv`.
- `config.yaml` is validated on load; a misspelt key fails at once.
- Messages go through `logging`; `-q` and `-v` set the level.

### Added
- Each stored mosaic keeps swissEO's cloud, terrain and snow classes, so the
  cloud rule is applied when reading and can change without a new download.
- Stable `parcel_id` in the parcel data (LV95 centroid and area).
- Per-field NDMI anomaly (`drought anomaly`) and swisstopo's VHI as a cross-check.
- Any GeoJSON polygon inside Switzerland as an AOI (`--aoi my_area.geojson`).
- Type checking (mypy) in CI and a pre-commit configuration.
