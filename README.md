# Reading a Drought from Space with Free & Open Satellite Data: Switzerland 2026

[![ci](https://github.com/floew2/drought-from-space-switzerland/actions/workflows/ci.yml/badge.svg)](https://github.com/floew2/drought-from-space-switzerland/actions/workflows/ci.yml)
![python](https://img.shields.io/badge/python-3.11%20%7C%203.13-blue)
![licence](https://img.shields.io/badge/licence-MIT-green)

**Built with** Python · GeoPandas · Shapely · rasterio · NumPy · pandas · SciPy ·
Matplotlib · Google Earth Engine Python API · pytest · GitHub Actions

**Background.** In summer 2026, Switzerland was hit hard by drought. At St. Gallen, rainfall to
mid-October reached only 61 % of the long-term median, the driest such period
in 31 years. Lake levels fell to record lows, the grass partly stopped regrowing, and
many farms had to feed their winter reserves in summer. This project looks at
the drought from space, using free and open satellite data
([Copernicus Sentinel-2](https://www.esa.int/Applications/Observing_the_Earth/Copernicus/Sentinel-2),
processed in [Google Earth Engine](https://earthengine.google.com/)) together
with open weather-station records. Its focus is grassland in three landscapes of
the canton of St. Gallen: the Linthebene, the St. Galler Rheintal (both valley
floors) and the Alpstein (a summering area). A secondary focus is making the
drought visible, field by field; the accompanying presentation also showed
browning across the Swiss Plateau and receding shorelines on Lake Constance,
which are not part of this repository.

![St. Galler Rheintal near Altstätten: Sentinel-2 true colour, 2019–2025 median against 2026, 15 June – 31 July](docs/figures/real/rheintal_s2_truecolour.jpg)

*St. Galler Rheintal near Altstätten, Sentinel-2 true colour. Left: median of
2019–2025, 15 June – 31 July. Right: 15 June – 31 July 2026. The
meadows browned field by field in summer 2026.*

**Data and method.** The 2026 season is compared with a 2019–2025 baseline,
using Sentinel-2 time series at 10 m (cloud-masked with Cloud Score+), the
cantonal agricultural parcel maps for 2017–2025 (meadows selected by their
official land-use (LNF) codes) and MeteoSwiss station data. Swiss meadows are
small: the median intensive meadow covers 0.67 ha and the median extensive one
just 0.19 ha, only about 67 and 19 pixels. Most pixels therefore touch a field
edge and mix grass with roads, hedges or crops. To keep that signal out, each
pixel's purity is modelled from the sensor's point-spread function (PSF): a
10 m detector box convolved with Gaussian optics, 12.7 m wide at half maximum,
derived from ESA's MTF specification. Only pixels that draw at least 90 % of
their signal from inside a meadow are used. That keeps about 80 % of the
intensive meadow area but only 45 % of the small, irregular extensive meadows;
the threshold trades monitored area against signal quality.

The analysis has three parts:

1. **Rainfall and water balance**, from MeteoSwiss Open Data.
2. **Grassland condition:** NDMI, a vegetation index sensitive to canopy water,
   is tracked over the pure pixels through each season and compared with the
   2019–2025 baseline.
3. **Mowing:** cuts are detected per field with a transparent rule set,
   validated against 143 hand-labelled field-years.

Details: [docs/METHOD.md](docs/METHOD.md).

![Three landscapes: Sentinel-2 imagery, field outlines, and pixel purity at three thresholds](docs/figures/real/mask_overview_intensive.png)

*Why pixel purity matters: intensive grassland at Linthebene, in the Rheintal
and at Alpstein. From left: Sentinel-2 false colour, the field
outlines, and purity at three thresholds (white pure, grey mixed, black
discarded). At purity ≥ 0 the pixels cover more than the field area, because
each pixel also sees past the field edge; at ≥ 0.9, 79–91 % of it remains.*

## Results

From May to September, the meadows' moisture index stayed 23–29 % below its
2019–2025 median in the Linthebene and the Rheintal, the lowest of the eight
years at both. For the Alpstein, a summering area, no season total is reported;
the maps below show a similarly widespread deficit there.

- Spring was normal. The meadows collapsed through July, bottomed out in
  mid-August and recovered from early September.
- St. Gallen recorded 690 mm of rain by 12 October against a median of 1,135 mm.
  The water balance P − ET0 stood at −74 mm against a median of +494 mm.
- From the hand labels, the last cut of the season came a median 12 days earlier
  in 2026 (15 of 22 fields). With 22 fields, that is not yet statistically
  significant (sign test p ≈ 0.13).
- **Mowing detection, validated like a model.** The rules were tuned on 42
  labelled field-years and then frozen; 78 field-years labelled later serve as a
  hold-out. On all 136 scorable field-years: precision 0.84, recall 0.90,
  F1 0.87 (0.90 without 2026). In 2026 precision drops to 0.64: 23 of 34 false
  detections fall in the late-August recovery, which mimics regrowth after a
  cut, as the synthetic data had predicted. A changepoint alternative scored
  F1 0.81 and was rejected.

![Precipitation at St. Gallen, grassland NDMI at Linthebene, and one meadow's mowing in 2025 against 2026](docs/figures/real/overview.png)

*Top left: cumulative rainfall at St. Gallen, 2026 against 1991–2020. Top right:
NDMI of intensive meadows at Linthebene, 2026 against 2019–2025. Bottom: one
meadow, two summers. Each sawtooth is a harvest; in 2026 the grass stops
regrowing after the July cut.*

![Per-field NDMI anomaly in 2026: St. Galler Rheintal, Linthebene and Alpstein / Churfirsten](docs/figures/real/field_anomaly_3sites.jpg)

*NDMI of each grassland parcel in 2026 against its own 2019–2025 median, in
standard deviations; red is drier than usual. From top: St. Galler Rheintal,
Linthebene, Alpstein / Churfirsten. Window May–August; at Alpstein, mid-July to
mid-September.*

**Next steps:** more hand labels for the mowing comparison, to test it with a
larger sample; validating the mowing detector outside the Linthebene; and
testing how the purity mask responds across ESA's MTF range.

More figures: [docs/gallery.md](docs/gallery.md). A step-by-step walkthrough:
[notebooks/walkthrough.ipynb](notebooks/walkthrough.ipynb).

## Quick start (synthetic data)

No Google Earth Engine account and no download of any data is needed to run the tools:

```bash
git clone https://github.com/floew2/drought-from-space-switzerland.git
cd drought-from-space-switzerland
make install        # pip install -e ".[dev]"   (uv: uv sync --extra dev)
make demo           # whole pipeline on synthetic data, ~90 s -> output/synthetic/
```

The synthetic inputs have the same formats as the real ones, so every step of
the pipeline runs on them. Figures go to `output/synthetic/figures/` and are
stamped SYNTHETIC; they show what each step does, not what actually happened in 2026.

## Running on real data

To run the tool on real satellite data, a Google Earth Engine account is needed
plus a set of vector parcel data from Canton St. Gallen over several years to
create the grassland mask.

| step                              | needs                                        | open? |
|-----------------------------------|----------------------------------------------|---|
| rainfall and water balance        | internet                                     | **yes**, MeteoSwiss Open Data |
| grassland mask                    | parcel polygons per year                     | **partly**: the latest year is free from [geodienste.ch](https://www.geodienste.ch/); earlier years are not open. With one year, every field is kept. |
| grassland condition               | the grassland mask + an Earth Engine project | **yes**, free for noncommercial use ([terms](https://earthengine.google.com/terms), [quotas](https://developers.google.com/earth-engine/guides/noncommercial_tiers)) |
| mowing scores (cutting frequency) | hand labels                                  | **no**; make your own with the labelling widget |

The parcel data used here came from the Canton of St. Gallen and is not
redistributed. Nothing in `data/` is committed.

```mermaid
flowchart LR
  P["parcel polygons"] --> M["mask build<br/>PSF purity, 10 m"]
  M --> U["upload + stamp-mask"]
  U --> E["Earth Engine exports"]
  S2["Sentinel-2 + Cloud Score+"] --> E
  E --> A["trajectories · sensitivity · cuts"]
  L["hand labels"] --> A
  W["MeteoSwiss"] --> R["precipitation"]
  A & R --> O["figures + results.json"]
```

<details>
<summary><b>Setup and commands</b></summary>

**Setup**

```bash
uv sync --extra dev --extra ee        # or: pip install -e ".[dev,ee]"
earthengine authenticate              # once; allow Drive access
```

Put your Earth Engine project in `config.local.yaml` (git-ignored):

```yaml
gee:
  project: <your-ee-project>
```

Exports go to your Google Drive and are downloaded from there: no Cloud Storage
bucket, no billing. Parcel shapefiles go to `data/nutzungen/<year>/`, hand
labels to `data/labels/`.

**Steps**

```bash
drought mask build                    # purity mask -> data/masks/grassland_mask_sg.tif
```

Upload the GeoTIFF once in the [Earth Engine Code Editor](https://code.earthengine.google.com/)
(Assets → New → GeoTIFF, name `grassland_mask_sg`, pyramiding **SAMPLE**), then:

```bash
drought ee stamp-mask                 # checks the asset against the local file
drought ee trajectory                 # site-mean index series
drought ee sensitivity                # purity-threshold sweep
drought ee cuts --labelled-only       # per-field series for the labelled fields
drought trajectory plot
drought sensitivity plot
drought precip --download
drought cuts score --site linthebene --labelled
drought cuts labels --site linthebene
drought results show                  # every headline number, refused if stale
```

Each `ee` command starts the export, waits, downloads and registers the file; an
interrupted one can simply be re-run. The trajectory export used about
25 EECU-hours of the free 150 a month. Every file carries a record of the
settings and mask that produced it, and loading refuses a mismatch.

To label mowing events for evaluation, in Jupyter:

```python
import drought as dr
df = dr.load_site("linthebene")
dr.labelling.label_widget(df, "linthebene", ev=dr.detect.detect_all(df), n=20)
```

All settings live in `config.yaml`; `drought --help` lists every command.

</details>

## Repository layout

```
drought-from-space-switzerland/
├── config.yaml                every setting
├── src/drought/               the analysis: mask, trajectory, precip, detect, evaluate, provenance
│   ├── engine/                Earth Engine exports, run from Python
│   └── synthetic/             synthetic stand-ins for every input
├── notebooks/
│   └── walkthrough.ipynb      the pipeline step by step
├── docs/
│   ├── METHOD.md              how and why
│   └── gallery.md             figures from the real data
└── tests/                     pytest; CI runs lint, tests, the demo and the notebook
```

## Data and licence

Code: MIT ([LICENSE](LICENSE)). Contains modified Copernicus Sentinel data
(2019–2026), processed in Google Earth Engine; Cloud Score+ (Pasquarella et al.
2023). Station data: MeteoSwiss Open Data (CC BY 4.0), *Quelle: MeteoSchweiz*.
Elevation: CGIAR-CSI SRTM 90 m. Field geometries: Nutzungsflächen, Kanton
St. Gallen, not redistributed.

## Author

Fabian Löw · [GitHub](https://github.com/floew2)
