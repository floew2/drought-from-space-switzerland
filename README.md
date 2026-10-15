# Reading a drought from space

Sentinel-2 and MeteoSwiss station data show that the 2026 summer drought in the
Canton of St. Gallen cut the season's grassland moisture signal by roughly a
quarter. Every site and management system ranked last of eight years.

![Linthebene, intensive meadows: 2026 NDMI against the 2019–2025 envelope](docs/figures/real/vi_trajectory_NDMI_intensive_linthebene.png)

*Linthebene, intensive meadows, site-mean NDMI. Real data, run of 14 October 2026.
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
| Season NDMI integral, intensive meadows | **−27 %** | **−23 %** |
| Season NDMI integral, extensive meadows | **−29 %** | **−28 %** |
| Rank among 2019–2026 | 1 of 8 | 1 of 8 |

The comparison is 2026 against the 2019–2025 median over May to September
(DOY 121–273), on the 103–114 days that every year observed. The deficit was not
season-long. Spring ran at or above normal; the index collapsed through July to
0.09–0.15 in mid-August and recovered sharply from early September.

The station data, an entirely independent source, agree. St. Gallen recorded
690 mm by 12 October, 61 % of the 1991–2020 median and the driest of 31 years.
The climatic water balance P − ET0 stood at −74 mm, against a median of +494 mm.

Every number above comes from the pipeline in this repository as it stands: the
box ⊛ Gaussian sensor response, no parcel size filter, the Earth Engine exports
run from Python, and every input checked against its provenance record. A
preliminary run with an older mask had given larger deficits at Linthebene (−33
and −35 %); [METHOD §3.4](docs/METHOD.md#34-result) says what changed and what
did not. The result does not depend on the residual-haze filter: with it off,
every series still ranks last, and three of four integrals move by under half a
point.

## What this repository is evidence of

The code is the smaller part. The judgement trail is in **[docs/METHOD.md](docs/METHOD.md)**;
in brief:

- **Which pixels to trust is derived, not assumed.** The Sentinel-2 response
  is modelled from the specified MTF: a detector box convolved with Gaussian
  optics. Pixel purity is the share of a cell's response that comes from
  inside a field. The threshold's cost was measured for each management
  system, and so was its effect on the answer. The 2026 anomaly grows steadily
  as the threshold tightens, by 3–4 % for intensive meadows and 9–15 % for
  extensive ones. The threshold in use, 0.90, is within 3.3 % of the strictest
  value tested, so mixed pixels make the reported drought signal slightly
  weaker, never stronger.
- **Invariants caught real bugs.** The observed area at purity > 0 must exceed
  the polygon area, because the response reaches past every boundary. Seeing
  it fall below 1 exposed a mask design that zeroed the halo, and the mask was
  rebuilt polygon-first. Matching parcels by centroid kept 43 %; matching by
  area overlap kept 86 % (81 % now that fields down to 0.05 ha are included).
- **The cut detector uses machine-learning methodology without a learned
  model.** It was tuned on 42 hand-labelled parcel-years and then frozen; the
  label set has since grown to 143 without the configuration changing. On the
  136 scorable parcel-years, F1 is 0.868, and 0.897 without 2026. In 2026,
  precision falls to 0.64: 23 of its 34 false positives fall between DOY 220
  and 260, on the drought recovery, which the synthetic data had predicted.
  So the last-cut comparison ends at DOY 228 for every year, and is also
  computed from the hand labels alone: the 2026 median is 12 days earlier,
  with 15 of 22 fields stopping sooner. That is an indication, not yet
  significant (sign test p ≈ 0.13). A changepoint alternative scored
  F1 0.805 on the same labels and was rejected. Both are kept in the code.
- **Numbers carry their provenance, and stale files are refused.** Every data
  file has a record of the settings and mask that produced it. Loading
  refuses a mismatch and names the field, which turns the project's worst
  recurring problem (results outliving the method that made them) into
  something enforced ([METHOD §6](docs/METHOD.md#6-bugs-found-by-checks)).
- **Packaging found things.** Packaging this repository turned up a
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
| `drought --synthetic cuts errors --site linthebene` | Where the detector errs: current year against the others, false positives by date. |
| `drought --synthetic cuts labels --site linthebene` | Last-cut shift from the hand labels alone, no detector (default end DOY 228). |
| `drought results show` | Every headline number from `results.json`, refused if any of it was computed under other settings or another mask. |
| `drought ee trajectory` / `sensitivity` / `cuts` | Runs the Earth Engine exports from Python: checks the mask asset is the local mask, waits, downloads, registers. |
| `drought gee render` | Writes `config.yaml`'s shared settings into the two interactive Code Editor scripts; CI fails if they drift. |
| `drought config show mask` | The effective configuration, after all overrides. |
| `make test` | One regression test per bug that actually happened, each also showing that the old logic fails it ([METHOD §6](docs/METHOD.md#6-bugs-found-by-checks)); plus the PSF derivation, CLI precedence and synthetic calibration. |

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
src/drought/engine/      the Earth Engine exports, run from Python (`drought ee …`)
gee/                     two interactive Code Editor scripts (map, hotspots)
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
