# Reading a drought from space

Sentinel-2 and MeteoSwiss station data show that the 2026 summer drought in the
Canton of St. Gallen cut the season's grassland moisture signal by a quarter to
a third. Every site and management system ranked last of eight years.

![Linthebene, intensive meadows: 2026 NDMI against the 2019–2025 envelope](docs/figures/real/vi_trajectory_NDMI_intensive_linthebene.png)

*Linthebene, intensive meadows, site-mean NDMI. Real data, run of {{RUN_DATE}}.
The input data cannot be published; the commands below run the same pipeline on
synthetic data.*

```bash
git clone https://github.com/floew2/drought-from-space-switzerland.git
cd drought-from-space-switzerland
make install        # pip install -e ".[dev]"
make demo           # whole pipeline on synthetic data, ~80 s -> output/synthetic/
```

---

## Results (real data)

| | Linthebene | St. Galler Rheintal |
|---|---|---|
| Season NDMI integral, intensive meadows | **−33 %** | **−23 %** |
| Season NDMI integral, extensive meadows | **−35 %** | **−30 %** |
| Rank among 2019–2026 | 1 of 8 | 1 of 8 |

The comparison is 2026 against the 2019–2025 median, integrated to DOY 245–250
over the 75–100 days that every year observed. The deficit was not season-long.
Spring ran at or above normal; the index collapsed through July to 0.09–0.15
in mid-August and recovered sharply from early September.

The station data, an entirely independent source, agree. St. Gallen recorded
501 mm by 10 August, 61 % of the 1991–2020 median and the driest of 31 years.
The climatic water balance P − ET0 ended at −84 mm, against a median of +324 mm.

All mask-based numbers above were produced by the mask as it was then built:
a single-Gaussian response model (σ 5.30 m) and a 0.25 ha parcel filter. The
filter has since been removed ([METHOD §2.1](docs/METHOD.md#21-which-fields-the-master-parcel-set)). The repository now defaults
to the model the method describes (see [METHOD §2.3](docs/METHOD.md#23-the-sensor-response)).
On synthetic parcels the model switch moves retained area by 1–3 points. The
real mask has not yet been rebuilt with either change.

## What this repository is evidence of

The code is the smaller part. The judgement trail is in **[docs/METHOD.md](docs/METHOD.md)**;
in brief:

- **Which pixels to trust is derived, not assumed.** The Sentinel-2 response
  is modelled from the specified MTF: a detector box convolved with Gaussian
  optics. Pixel purity is the share of a cell's response that comes from
  inside a field. The threshold's cost was measured for each management
  system.
- **Invariants caught real bugs.** The observed area at purity > 0 must exceed
  the polygon area, because the response reaches past every boundary. Seeing
  it fall below 1 exposed a mask design that zeroed the halo, and the mask was
  rebuilt polygon-first. Matching parcels by centroid kept 43 %; matching by
  area overlap kept 86 %.
- **The cut detector uses machine-learning methodology without a learned
  model.** It was tuned on 42 hand-labelled parcel-years and then frozen.
  Validated on all 65, F1 went from 0.911 to 0.912, and recall per year was
  0.81–1.00. A changepoint alternative was implemented, scored on the same
  labels, and rejected at F1 0.871. Both are kept in the code.
- **Numbers carry their provenance.** Packaging this repository turned up a
  mismatch between the response model on the slides and the one in the code,
  and a filter bug in the parcel selection. Both are documented and measured,
  and neither was fixed silently ([METHOD §6](docs/METHOD.md#6-bugs-found-by-checks)).

## Running it

Everything runs through one command, `drought` (see `drought --help`):

| | |
|---|---|
| `drought demo` | Mask build, area accounting, trajectories, precipitation, cut detection and scoring, all on synthetic data. Exits non-zero if a calibrated value misses its target or a physical invariant fails. `make demo` runs the same. |
| `drought --synthetic mask check` | Area cost of each purity threshold, and the PSF-model comparison. |
| `drought --synthetic cuts score --site linthebene` | The frozen detector and the changepoint alternative, scored against the labels. |
| `drought config show mask` | The effective configuration, after all overrides. |
| `make test` | Unit tests: PSF derivation and the response width the code applies, tiled = single-pass purity, the master-set funnel, CLI precedence, synthetic calibration. |

Settings come from `config.yaml`. Precedence, lowest to highest:
`config.yaml` < `--config FILE` < `--synthetic` < `--set KEY=VALUE`. For
example, `drought --synthetic --set mask.psf_model=gauss mask build` rebuilds
the mask with the old response model. Only existing keys are accepted. An
override that would arrive after a module has already read the setting is
refused rather than silently ignored.

The synthetic inputs are generated deterministically by `drought.synthetic`.
They match the real data's schemas, its parcel-size and shape distributions,
its cadastral redraws, splits and merges, and its headline integrals. The demo
prints every synthetic number beside the real one, marked *calibrated* or
*emergent*. A calibrated number reproducing its target shows the pipeline is
wired correctly. It is not new evidence about the drought. Every
synthetic figure is stamped SYNTHETIC.

To run on the real data, see [docs/WORKFLOW.md](docs/WORKFLOW.md).

## Repository map

```
config.yaml              every setting, one place
src/drought/             the analysis (mask, trajectory, precip, detect, evaluate, ...)
src/drought/synthetic/   calibrated stand-ins for every input; targets.py lists the real numbers
gee/                     Earth Engine extraction scripts
notebooks/analysis.ipynb the working notebook (outputs stripped)
docs/METHOD.md           how and why — start here
docs/gallery.md          rendered figures: real results, model, synthetic demo
docs/WORKFLOW.md         runbook for the real data
tests/                   pytest; CI runs lint, tests and the demo
```

## Data and licence

Code: MIT ([LICENSE](LICENSE)). Contains modified Copernicus Sentinel data
(2019–2026), processed in Google Earth Engine. Station data: Source: MeteoSwiss
(CC BY 4.0). Field geometries and land-use codes: Nutzungsflächen, Kanton
St. Gallen. These are **not redistributed** in any form.
