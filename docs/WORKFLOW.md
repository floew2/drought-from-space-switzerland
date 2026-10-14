# Workflow — running the analysis on the real data

This is the runbook for the real inputs, which are not in this repository. For
what the analysis found and why it is built this way, see [METHOD.md](METHOD.md);
to run it on synthetic data, `make demo` (see the [README](../README.md)).

---

<!-- previous title: Drought 2026 — Swiss grassland and water -->

Sentinel-2, MODIS and MeteoSwiss analysis of the 2026 Swiss drought: grassland
condition, mowing behaviour, and lake water extent.

The work splits into two kinds of thing, and it matters which is which:

| | runs where | produces | goes into |
|---|---|---|---|
| **Visuals** | GEE Code Editor | map screenshots, GeoTIFFs | slides |
| **Analysis** | GEE → CSV → notebook | figures with numbers | slides + write-up |

Some Earth Engine scripts only make pictures. Others extract data you then
download and analyse locally. Mixing them up wastes compute, and compute here is
not free: **~40,000 EECU-seconds per site-year** for the per-parcel extraction.

---

## 0. Setup

```bash
make install        # pip install -e ".[dev]"   (with uv: uv sync --extra dev)
drought demo        # the whole pipeline on synthetic data -> output/synthetic/
drought --help      # every command
```

For the real data, the same commands run without `--synthetic`, against the
paths in `config.yaml`:

```bash
drought mask build                    # needs data/nutzungen/
drought mask check
drought trajectory plot --site linthebene
drought precip --download             # or reads data/station/ if cached
drought cuts score --site linthebene  # needs data/gee/<site>/ and data/labels/
drought run-all                       # every step that has its input; reports the rest
```

The real inputs (cantonal parcel data and everything derived from it) are not
in this repository and cannot be published. `make demo` runs the same code on
synthetic data calibrated to the real distributions; see
`src/drought/synthetic/`. For the real analysis, put the data in `data/` (see
`data/README.md`) and use the notebook:

```bash
pip install -e ".[notebook]"
jupyter lab notebooks/analysis.ipynb
```

All settings live in `config.yaml`. Change them **there**, never in the
notebook, so any figure can be traced back to what produced it.

```
config.yaml            all settings
notebooks/analysis.ipynb
src/drought/           the code
src/drought/engine/    the Earth Engine exports, run from Python
gee/                   two interactive Earth Engine scripts (maps, screenshots)
data/nutzungen/        yearly Nutzungsflaechen shapefiles  <- you provide
data/gee/<site>/       Earth Engine exports                 <- `drought ee` downloads
data/labels/           hand-labelled reference data
data/masks/            grassland mask GeoTIFF
output/figures/<site>/
output/tables/         provenance JSON beside every result
```

---

## 1. Earth Engine — run from here, not from the browser

The three data exports run from Python (`drought.engine`, earthengine-api).
Every setting comes from `config.yaml`; nothing is pasted or edited by hand.

**Once:**

```bash
uv sync --extra dev --extra ee
earthengine authenticate        # allow Drive access when asked
```

Put the project in `config.local.yaml`, which git ignores and which may set
nothing but the project and the Drive folder:

```yaml
gee:
  project: <your-ee-project>
```

No Cloud Storage bucket and no billing: exports are batch tasks to Google
Drive (folder `gee.drive_folder`), downloaded from there with the Earth Engine
login. If Drive refuses the download, enable the Google Drive API in the
project (Cloud Console → APIs & Services; free).

**Every run:**

```bash
drought mask build              # the mask and its provenance record
```

Upload `data/masks/grassland_mask_sg.tif` in the
[Code Editor](https://code.earthengine.google.com/): Assets → NEW → GeoTIFF,
asset name `grassland_mask_sg`, pyramiding policy **SAMPLE** (the default,
MEAN, averages purities at coarser scales). Replace the old asset. Then:

```bash
drought ee stamp-mask           # asset == local mask? then stamp it with the mask's hash
drought ee check                # project, Drive, and the stamp
drought ee trajectory --smoke   # one site, two years: read the EECU cost first
drought ee trajectory           # -> data/gee/grassland_trajectory_multi_p90_cfg<hash>.csv
drought ee sensitivity          # -> data/gee/purity_sensitivity_<site>_<system>_cfg<hash>.csv
drought ee cuts --smoke         # current year, 100 parcels: read the EECU cost first
drought ee cuts --labelled-only # only the hand-labelled parcel-years -> data/gee/<site>/labelled/
drought cuts score --site linthebene --labelled   # re-score the frozen detector on them
drought ee cuts                 # -> data/gee/<site>/cuts_<site>_<year>_cfg<hash>.csv
```

Each command starts the tasks, waits (printing state changes and EECU cost),
downloads the CSVs from Drive and registers them. A command interrupted while
waiting can simply be run again: finished exports are downloaded, not repeated.

**Why nothing stale gets through:**

- **The asset must be the local mask.** `stamp-mask` reads 400 cells of the
  uploaded asset (fixed seed, mostly grassland) and compares all four purity
  values with the local GeoTIFF. Only if every one agrees to 1e-6 does it store
  the mask's content hash on the asset (`drought_mask_sha256`). A wrong file,
  a stale upload, a shifted grid or swapped bands all fail. Every export then
  refuses to start unless the stamp equals the local mask.
- **Exports are named with a hash** of every Earth Engine setting plus the mask
  (`_cfg<hash>`). Registration refuses a name from other settings and names the
  setting that changed.
- **Registration checks the content**: the management categories, the years,
  and no geometry column in a cut export. Every analysis step then refuses an
  unregistered or stale file (`provenance.require`).

The smoke runs land in `output/ee_smoke/` and are not registered.

**Cost.** The per-parcel cut extraction is by far the most expensive export.
The smoke run (2026, 100 parcels) cost 76,000 EECU-s at Linthebene and
18,000 at Rheintal. Linthebene lies on two MGRS tiles and has the most
acquisitions. The full sample (250 parcels, eight years, two sites) is
therefore of the order of a million EECU-s, several times the noncommercial
Community tier's 150 EECU-hours a month. `--labelled-only` re-scores the
detector on the labelled parcel-years alone. Cleaning works per parcel, so
their series are the same as inside the full sample.

### For screenshots only — nothing to download
### For screenshots only — nothing to download

| script | what | how to use |
|---|---|---|
| `gee_switzerland_rgb.js` | National MODIS true-colour pair, 2026 vs a ten-year baseline, 500 m | Set `YEAR_DRY` / `BASELINE_YEARS`, run, screenshot the two layers. Also exports GeoTIFFs to Drive for print resolution. |
| `gee_s2_hotspots.js` | Sentinel-2 zoom-ins at 10 m: Untersee, Horn, Rheintal, Linthebene. RGB pairs plus lake water extent by two methods | Set `ACTIVE_SITE`, run, screenshot. Run **once per site**. |
| `gee_grassland_map.js` | Per-parcel grassland condition anomaly (NDMI z-score), choropleth | Set `ACTIVE_SITE`, run, screenshot. Needs the mask asset from step 2. |

These run in the [Code Editor](https://code.earthengine.google.com/). Their
shared settings are generated from `config.yaml`, and the committed copies have
a placeholder project, so paste from a rendered copy:
`drought --set gee.project=<project> gee render --out output/gee`.

Stretches are **pinned** — identical `min`/`max` on baseline and current year,
no per-scene levelling. That is what makes the two panels comparable, so don't
"improve" the contrast on one of them.

---

## 2. Build the grassland mask (once, then rarely)

Everything about grassland depends on knowing which pixels are grassland, and on
that being true for **every** year, not just the current one.

1. Put the yearly Nutzungsflaechen shapefiles in `data/nutzungen/`, 2017–2025,
   named `nutzungen_<year>.shp` with their `.shx`, `.dbf` and `.prj`. One folder
   per year is fine; the search is recursive. A year without a `.prj` is refused.
2. Build:
   ```bash
   drought mask build      # prints the parcel funnel; writes the mask + its provenance record
   drought mask figures    # parcel_size, area_check, purity explainer
   ```
   GDAL warns about a handful of self-touching rings while reading; the build
   repairs them (`make_valid`) and prints how many per year.
3. Upload the resulting `data/masks/grassland_mask_sg.tif` to GEE as an **image
   asset** named `grassland_mask_sg`.

The GEE scripts read that asset. They never see the polygons — which is the
point: rasterising parcels inside a per-image loop once burned 137,000
EECU-seconds before it was killed.

**Why persistence matters.** A mask built from one year says "grassland in
2025", and the analysis then reads a 2019 reflectance value through it. A parcel
that was maize in 2019 would put an arable curve into the grassland baseline.
The mask therefore keeps only land that was grassland in **every** year, and
exports both variants so the choice can be tested rather than assumed.

---

## 3. Workflow, in order

### A. Precipitation and water balance — no GEE

Notebook chapter 2:

```python
dr.precip.main()
```

Downloads MeteoSwiss station data directly. Current year and the 1991–2020
reference come from the **same series** — mixing a reanalysis current year with
published station normals would put the model's bias straight into the
"deficit".

Two panels: cumulative precipitation, and cumulative P − ET0 (Hargreaves). The
second is an index of atmospheric demand against supply, **not** a soil water
balance.

### B. Imagery for slides — GEE only

Run `gee_switzerland_rgb.js`, then `gee_s2_hotspots.js` once per site.
Screenshot. Nothing to download unless you want print-resolution GeoTIFFs.

**Read the printed solar zenith angles** in the MODIS script. Terra's overpass
has drifted since 2020 (28.5° → 41.1° over the compositing window), so part of
the darkening in the national panels is lower illumination, not drought. That
belongs in the caption. Sentinel-2 holds a stable overpass, which is why every
impact figure uses Sentinel-2 rather than MODIS.

### C. Grassland condition — trajectory

1. `drought ee trajectory` (§1).
2. `drought trajectory plot --site linthebene`, or notebook chapter 4:
   ```python
   dr.trajectory.plot(site="linthebene", index="NDMI")
   ```

The trajectory panel is the trustworthy one — it shows *when* the divergence
happened and depends on no interpolation choices. The season integral is
sensitive to how gaps are handled, which is why it is the second panel.

### D. Grassland condition — map

Needs the mask asset. Run `gee_grassland_map.js`, screenshot.

### E. Cut detection

1. `drought ee cuts --smoke`, read the EECU cost, then `drought ee cuts`.
   Cost predictions here have been wrong more than once. The tasks run newest
   year first, so a run cut short still holds 2026.
2. Notebook chapter 5:
   ```python
   df = dr.load_site("linthebene")
   ev = dr.detect.detect_all(df)
   dr.figures.figure_last_cut(df, "linthebene")
   ```

**Which metric to trust**, in order:

- **Last cut of the season** — the detector needs to find *one* cut per
  parcel-year, not all of them, and the bias runs the right way: the drought
  year has the most images, so a late cut is easier to detect there.
- **Late-season share** — a ratio, so sparse sampling loses cuts from both
  halves roughly equally.
- **Cut count** — *not comparable between years* without density matching.
  Sensitivity rises from ~1.3 cuts at <15 clear images to ~2.9 at 35–45, and
  drought summers are sunny.

### F. Reference data and validation — do not skip

The detector is a heuristic. Without hand labels you cannot tell whether it
works, and it has been confidently wrong before.

1. Label ~20 parcel-years per site (chapter 6.1):
   ```python
   dr.labelling.label_widget(df, SITE, ev=ev, n=20, seed=0)
   ```
   Click to add at the nearest observation, click a marker to remove. **H**
   hides the detector's proposals — label blind first, then reveal and compare.
   Labels copied from the detector measure nothing.

   Place each label at the **trough**: the first point clearly down.

2. Save, then score:
   ```python
   labels = dr.labels.check_labels(labels, df)
   dr.labels.save_labels(labels, SITE)
   dr.evaluate.score_by_year(labels, df)
   ```

**Read `score_by_year`, not the overall F1.** The headline compares years, so
what threatens it is not accuracy but accuracy that *differs* by year.

Labels persist in `data/labels/cut_labels_<site>.json` and load automatically
next session, with already-labelled parcel-years excluded from new samples.
**One file per site** — a detector validated on Linthebene has not been
validated on Rheintal.

---

## 4. Module reference

| module | what |
|---|---|
| `config` | loads `config.yaml`, resolves site-aware paths, `provenance()` |
| `io` | reads GEE exports (both wide and long generations) |
| `clean` | per-parcel coverage filter, residual-haze filter |
| `detect` | the heuristic cut detector (frozen configuration) |
| `changepoint` | alternative method — tested, rejected, kept as the record |
| `metrics` | late share, density-equalised counts, last-cut, paired comparisons |
| `evaluate` | scoring against hand labels, per-year, tuning, `assert_wired` |
| `labels` | save / load / validate reference data |
| `labelling` | interactive click-to-label widget (pure HTML/JS, no deps) |
| `figures` | all plots, written to `output/figures/<site>/` |
| `precip` | MeteoSwiss cumulative precipitation and P − ET0 |
| `trajectory` | site-mean VI trajectory and season integral |
| `mask` | builds the grassland mask from yearly shapefiles |
| `synthetic` | calibrated stand-ins for every input, for `make demo` |
| `demo` | the whole pipeline on synthetic data, self-checking |
| `gallery` | refreshes the committed figures in `docs/figures/` |

---

## 5. Data sources

Copernicus Sentinel-2 MSI Level-2A, harmonized, 10 m — contains modified
Copernicus Sentinel data 2019–2026. NASA LP DAAC MODIS/Terra MOD09A1 v061.
Cloud masking: Cloud Score+ (Pasquarella et al. 2023, CVPRW,
doi:10.1109/CVPRW59228.2023.00206). Land cover: MCD12Q1 v061, ESA WorldCover
v200, Google Dynamic World V1. Elevation: CGIAR-CSI SRTM 90 m v4. Boundaries:
FAO GAUL (2015). Parcels: Nutzungsflächen, Kanton St. Gallen, annual
declarations 2017–2025. Precipitation: MeteoSchweiz Open Data, SwissMetNet
station St. Gallen, daily values (`rre150d0`), CC-BY — *Quelle: MeteoSchweiz*.
Projection: CH1903+ / LV95 (EPSG:2056).
