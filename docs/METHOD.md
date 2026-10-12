# Method — how the analysis works, and why it is built this way

This document records how the analysis was built and why. The interesting part
is the decisions: what was measured before anything was chosen, which checks
caught which errors, and what was rejected. The code says *what*; this file says
*why*, and where the evidence is.

All real-data numbers are from runs on cantonal data that is not in this
repository. Where a number comes from an earlier version of the pipeline, it is
labelled as such. Synthetic numbers are always labelled synthetic.

**Contents**
1. [Question and data](#1-question-and-data)
2. [The grassland mask](#2-the-grassland-mask)
3. [Grassland condition](#3-grassland-condition)
4. [Precipitation and water balance](#4-precipitation-and-water-balance)
5. [Cut detection](#5-cut-detection)
6. [Bugs found by checks](#6-bugs-found-by-checks)
7. [The synthetic data, and what it can and cannot show](#7-the-synthetic-data-and-what-it-can-and-cannot-show)
8. [Open items](#8-open-items)

---

## 1. Question and data

How hard did the 2026 summer drought hit grassland in the Canton of St. Gallen?
Was the effect visible in farming behaviour (mowing) as well as in vegetation
condition? Three sites cover the range: two valley-floor grassland areas
(Linthebene, St. Galler Rheintal) and one alpine summering area
(Alpstein / Churfirsten).

| input | what | role |
|---|---|---|
| Sentinel-2 L2A, harmonized, 10 m, 2019–2026 | NDVI, kNDVI, NDMI, NDRE; Cloud Score+ masking | condition and cuts |
| Nutzungsflächen 2017–2025 | annual declarations of field outlines and land-use codes | which pixels are grassland, and of which kind |
| MeteoSwiss SwissMetNet, St. Gallen, 1991–2026 | daily precipitation and temperature | independent check |

The baseline is 2019–2025. That baseline includes the 2022 drought summer, so a
2026 anomaly measured against it is conservative: a floor, not a ceiling.
Everything is in CH1903+ / LV95 (EPSG:2056).

---

## 2. The grassland mask

A Sentinel-2 pixel over Swiss grassland is often not purely grassland, and when
it is, it may not be the same kind of grassland. The mask answers one question
per pixel: *how much of what the sensor records here comes from inside one
management category of grassland?*

### 2.1 Which fields: the master parcel set

The mask is built **polygon-first**. A fixed, countable set of fields is
selected, and purity is computed from those polygons alone.

- A field qualifies if it was grassland in every year from 2017 to 2025 **and**
  held the same management category in every year. A field that was extensive
  until 2024 and intensive since would otherwise carry an extensive history into
  an intensive baseline.
- **Parcels are matched between years by area overlap, not by identity.** A
  2025 field "held category C in year Y" if at least 80 % of its area falls
  inside year Y's grassland of category C. Outlines are redrawn, split and
  merged between annual deliveries, so every method that tried to follow an
  individual parcel failed. Centroid matching kept 43 %; overlap matching keeps
  86 %.
- **Code 930 (Sömmerung) is exempt** from the persistence test, rather than
  failing it. It appears only in the most recent delivery, so there is nothing
  to check it against. That is a schema change, not a change of land use.
- **Persistence is a property of the parcel set and is never multiplied into
  the raster.** [§6](#6-bugs-found-by-checks), bug 1, explains why this matters.
- **There is no size filter.** Whether a field is large enough for Sentinel-2
  is decided by pixel purity (§2.3–2.4), and by nothing else. A field that
  yields no cell above the purity threshold drops out on its own. The only floor
  is a **sliver floor at 0.05 ha**, applied in every year. It removes cadastral
  offcuts (digitising leftovers and boundary slivers) that would corrupt the
  area-overlap matching. It is named for that job, not for observability.

An earlier version also applied a 0.25 ha filter, in every year. That was a
second mechanism doing purity's job, and running it on past deliveries failed
stable fields ([§6](#6-bugs-found-by-checks), #11). Size is also a poor proxy
for what the sensor can resolve. Here is the share of isolated rectangular
fields, at random orientation, that yield at least one cell with purity ≥ 0.9
under the box ⊛ Gaussian response:

| shape index | 0.05 ha | 0.10 | 0.15 | 0.20 | 0.25 | 0.30 | 0.40 | 0.50 |
|---|---|---|---|---|---|---|---|---|
| 4.0 (square) | 40 % | 100 % | 100 % | 100 % | 100 % | 100 % | 100 % | 100 % |
| 5.72 (median intensive) | 0 % | 0 % | 0 % | 100 % | 100 % | 100 % | 100 % | 100 % |
| 6.29 (median extensive) | 0 % | 0 % | 0 % | 5 % | 95 % | 100 % | 100 % | 100 % |
| 8.0 | 0 % | 0 % | 0 % | 0 % | 0 % | 0 % | 72 % | 100 % |

A compact 0.1 ha field is observable; a 0.3 ha strip with shape index 8 is
not. A fixed area cut gets both wrong. Purity gets both right, by definition.

Real funnel, full extent. **This was measured under the old 0.25 ha filter, so
it is provisional.** It will be re-measured on the next real run.

```
 40,205  grassland parcels >= 0.25 ha in 2025
  3,592  dropped — not grassland in enough years
  2,027  dropped — changed management category
  2,137  exempt (code 930, latest delivery only)
 34,586  remain (86.0 %), 80,850 ha
```

The analysis keeps two category systems, and they do not nest:

| grouping | LNF codes | the question it answers |
|---|---|---|
| intensive / extensive | 613, 650 / 611, 612, 634, 635, 693, 694 | the farming system: 4–6 cuts from early May, against about 2 cuts with none before 15 June |
| mown / grazed | from `config.yaml` | whether a cut signal exists at all |

Intensive and extensive meadows are never pooled. Averaging two cutting regimes
gives a curve that belongs to neither.

### 2.2 Why purity matters: the fields are small and mostly edge

Real master set, measured under the old 0.25 ha filter (provisional):

| category | parcels | median ha | < 0.5 ha | median shape index | 10 m pixels in a median field |
|---|---|---|---|---|---|
| intensive | 23,758 | 1.11 | 21 % | 5.72 | 111 |
| extensive | 4,348 | 0.45 | 58 % | 6.29 | 45 |
| mown | 28,112 | 0.93 | 26 % | 5.80 | 93 |
| grazed | 6,474 | 0.99 | 25 % | 6.61 | 99 |

The shape index is perimeter / √area, where 4.0 is a perfect square. A median
intensive field holds about a hundred 10 m pixels, which sounds comfortable.
But a shape index near 6 means a long, narrow field, and much of it lies within
a pixel or two of its boundary.

This is also why the grid is **10 m**, the Sentinel-2 native grid, rather than
the 20 m used earlier. A 0.44 ha field covers 11 cells at 20 m and 44 at 10 m,
so edge effects are resolved rather than assumed.

### 2.3 The sensor response

A pixel does not see a sharp 10 × 10 m box. The system point-spread function
also draws signal from beyond the nominal footprint. Purity is therefore defined
on the *response*, not on area. Parcels are rasterised at 2.5 m and convolved
with the modelled response, and each 10 m cell gets the share of its response
that comes from inside the category's parcels.

The response is derived from the system modulation transfer function (MTF) at
Nyquist:

| quantity | value | how |
|---|---|---|
| system MTF at Nyquist | 0.25 | ESA MSI specification (citation pending, [§8](#8-open-items)) |
| detector MTF at Nyquist | 0.637 | the detector footprint is a box one pixel wide: sinc(½) |
| optics MTF at Nyquist | 0.393 | 0.25 / 0.637 |
| σ, Gaussian optics | 4.35 m | from MTF(f) = exp(−2π²σ²f²), f = 1/(2·10 m) |
| FWHM of the response | 12.7 m | box ⊛ Gaussian, measured numerically |

A Gaussian with σ 4.35 m on its own would be only 10.2 m wide. The 12.7 m exists
only after convolution with the detector box. The two must never be stated as
one Gaussian.

**What built the published numbers.** The code that built the mask did not
implement this model. `purity_psf()` used a single Gaussian fitted to the whole
system MTF, with σ 5.30 m. It read each cell as the mean of the 2 × 2 sub-cells
around its centre. The deck and `figure_psf()`, meanwhile, showed box ⊛
Gaussian. Every real number in §2.4 comes from the single-Gaussian mask.

Both models reach MTF 0.25 at Nyquist, and their second moments agree to 1.5 %
(5.30 m against 5.22 m). `mask.measured_response_fwhm()` passes a line source
through the code and measures the width it actually applies. Under that
measurement the old implementation's 2 × 2 read-out adds about 5 m of its own
span, so its effective response is about 12.95 m, not the nominal 12.48 m.

The repository now defaults to box ⊛ Gaussian (`mask.psf_model: box_gauss`),
with the box applied as a full-cell average. The old model is still selectable,
and `mask.psf_model_check()` builds both on the same parcels. On the synthetic
landscape:

| category | area ratio at ≥ 0.9, box ⊛ G | single G | difference |
|---|---|---|---|
| intensive | 0.783 | 0.769 | +1.4 points |
| extensive | 0.482 | 0.450 | +3.2 points |
| mown | 0.777 | 0.764 | +1.4 points |
| grazed | 0.890 | 0.882 | +0.7 points |

So the model is not a lever on the result, but it is not nothing either. The
real mask has not been rebuilt, and the published numbers carry the old model's
label.

**A trap when switching.** Swapping σ from 5.30 m to 4.35 m while keeping the
centre read-out gives Gaussian(4.35) sampled at the centre. That response is
about 11 m wide, narrower than either real model, so it under-weights the halo.
A test (`test_implemented_response_has_the_stated_width`) measures the response
the code actually applies and fails on exactly this mistake.

Approximations to state with the result: the model is nadir-only and isotropic;
real optics are Airy-like rather than Gaussian; the real response is slightly
wider along-track because of motion during integration; and ESA publishes MTF
per band, while one representative value is used here.

### 2.4 What the threshold costs

Real data, master set clipped to the raster window. These numbers were measured
with the single-Gaussian model and the old 0.25 ha filter, so they are
provisional:

| category | parcels | polygon ha | observed / polygon area at purity > 0 | at ≥ 0.9 |
|---|---|---|---|---|
| mown | 23,662 | 40,116 | 1.54 | 0.84 |
| grazed | 5,268 | 18,418 | 1.49 | 0.88 |
| intensive | 20,094 | 37,996 | 1.54 | 0.84 |
| extensive | 3,563 | 2,112 | 3.40 | 0.56 |

At the threshold of 0.9, intensive meadows keep 84 % of their mapped area. The
small, convoluted extensive meadows keep 56 %, which is the price of measuring
them without their neighbours mixed in.

**The ratio at > 0 must exceed 1.** The response reaches past every boundary,
so the cells it touches always cover more ground than the polygons. This is the
check that caught the worst bug in the project ([§6](#6-bugs-found-by-checks),
bug 1), and the demo enforces it. It is a check, not a measurement. Its exact
value depends on where the numerical tail of the response is cut off: rounding
purity to six decimals moved the synthetic intensive ratio from 1.44 to 1.41.
Crossing 1 is not a target either; that happens where roughly half the retained
cells are contaminated.

### 2.5 Is the threshold a lever on the answer?

The purity threshold is applied in Earth Engine, not baked into the raster, so
it can be swept at the cost of one extraction. On the **previous 20 m mask**,
the 2026 NDMI anomaly moved by under 1 % of its own size across the whole range
0.05–0.95, at two independent sites. If that holds at 10 m, the threshold
changes how many pixels are measured but not what they say. **It has not yet
been re-run on the 10 m mask** ([§8](#8-open-items)).

That earlier sweep had a design flaw. It admitted only parcels with at least
one cell of purity 0.90, so composition stayed fixed while the threshold
moved. That isolates contamination cleanly, but only for large, compact fields,
which are the fields contamination barely touches. The re-run admits every
parcel (`gee.strict_purity: 0.0`). It answers the question the headline
depends on: would the reported anomaly change under another threshold, for the
grassland it is reported for? At high thresholds small fields drop out, so a
drift there could be contamination or composition. If the curve drifts, the
0.90 run separates the two. At Alpstein the sweep uses the summering window,
DOY 196–258, rather than May–August. The sweep uses the same Cloud Score+
threshold as the trajectory script (0.60), so it tests the pixels behind the
reported number.

---

## 3. Grassland condition

### 3.1 Site-mean trajectories, and what is thrown away

For each site, management system, year and acquisition date, Earth Engine
returns the mean index over mask pixels with purity ≥ 0.9, plus the pixel count.
Before any curve is drawn, three filters run:

- **Tile duplicates.** Two Sentinel-2 tiles see the same ground on the same day.
  The index is averaged, weighted by pixel count. The count itself takes the
  maximum, not the sum: summing would double-count the overlap and make every
  ordinary date look half-covered.
- **Partial coverage.** A date is dropped if less than 80 % of the site's mask is
  clear. A partly clouded scene samples a *different set of fields*, not a random
  subsample. It was the leading cause of large day-to-day jumps.
- **Cloud Score+ threshold.** A pixel is used only if its Cloud Score+
  `cs_cdf` (a per-pixel clear-sky score, 0–1) is at least 0.60, Google's
  suggested default within its 0.50–0.65 range. The cut extraction is the one
  exception, at 0.50 (§5.2). The two values were inherited rather than chosen;
  each is now set once in `config.yaml` and rendered into the scripts.
- **Residual haze.** Cloud Score+ lets thin haze through, and haze leaves a
  fingerprint ordered by wavelength. Scattering inflates the short-wavelength
  bands, so red is hit hard, red-edge less, and SWIR hardly at all. One
  Rheintal date in 2026 lost 0.197 NDVI, 0.086 NDRE and 0.028 NDMI over two days
  at full coverage, which no vegetation change does. Dates are flagged where NDVI
  falls well below its local level while NDMI does not follow (robust z-scores,
  so no absolute threshold has to be guessed). About 17 % of dates go, spread
  evenly across years, which shifts the 2026 anomaly slightly further negative.

The current year keeps its raw observations as dots beside the smoothed curve,
so smoothing never hides what was measured.

### 3.2 The season integral

The integral of NDMI over a fixed window is an established proxy for season
production, and the window is DOY 121–273 (1 May – 30 September). Two earlier
versions of this calculation were wrong, and each was caught by a
contradiction in its own output:

- **A data-driven window.** The window was first the intersection of the dates
  every year covered. One cloudy September in a single baseline year closed the
  window at mid-August for every year, which is exactly where the 2026 deficit
  begins. It still accumulated two months of an unusually strong spring. The
  result was *positive* integrals alongside *negative* growing-season means.
  The window is now fixed.
- **Cloud luck.** A gap counted as zero growth, so a well-observed year beat a
  poorly observed one regardless of vegetation. The symptom was a curve 15 %
  below the median that ranked only 4th of 8. Every year is now integrated over
  exactly the days that **all** years observed: 75–100 of the 153 in the window,
  depending on site and management system.

The integral stops at the current year's last observation, so it is a
season-to-date comparison, not a full season. Gaps longer than 30 days are left
empty rather than interpolated across. A three-week hole would otherwise
invent exactly the kind of shape a sceptical reader would go looking for.

### 3.3 Why NDMI, read as an anomaly

NDMI (NIR against SWIR) responds to canopy water and keeps responding over dense
swards. NDVI saturates near 0.9 there. NDMI has no accepted absolute scale,
because it mixes water content with canopy density: the same leaf water reads
anywhere from 0.10 to 0.39 depending on leaf area. So it is only ever read
against each site's own history, never as an absolute value. All four indices
were computed. If they agree on sign and rough size, the result does not depend
on which one is shown.

### 3.4 Result

| | Linthebene | St. Galler Rheintal |
|---|---|---|
| intensive | −33 % | −23 % |
| extensive | −35 % | −30 % |

Every one of the four was the lowest of eight years. Spring ran at or above
normal; the index collapsed through July to 0.09–0.15 in mid-August and
recovered sharply from early September. That is a severe midsummer drought that
broke, not a failed year. The two management systems have very different
cutting regimes, so their agreement rules out mowing as the explanation.
Rheintal is consistently milder, as its position near the Rhine and its
shallower groundwater would predict.

---

## 4. Precipitation and water balance

The current year and the 1991–2020 reference come from the **same** station
series (MeteoSwiss SwissMetNet, St. Gallen, daily). Using one series means a
difference cannot be an instrument or source change. Days are aligned by
calendar date, with 29 February dropped, so the current year is cut at exactly
the same day as every baseline year. A monthly cumulative sum would have given a
partial-month artefact. The figure shows the envelope of the 30 individual
baseline years, so a reader can see whether 2026 lies outside the historical
range, not merely below average.

To 10 August 2026: **501 mm against a median of 826 mm**. That is 61 % of the
median and the driest of 31 years in the comparison window.

Precipitation alone understates a hot drought, so the second panel is the
climatic water balance P − ET0, with ET0 from Hargreaves. This is an index of
atmospheric demand against supply, **not** a soil water balance. It ended at
**−84 mm against a median of +324 mm**. Because the quantity crosses zero, it is
reported as a difference, never as a percentage. "−26 % of normal" would be a
serious misreading of a 408 mm shortfall.

---

## 5. Cut detection

### 5.1 Why per field, and why NDMI

A site mean averages hundreds of fields cut on different dates, which smooths
the mowing sawtooth away almost entirely. A cut is a field-level event, so it is
measured per field. Each field is buffered inward by 15 m, and only fields of at
least 0.5 ha are sampled.

The detector reads NDMI and uses NDVI to confirm. Over dense grassland NDVI falls
only about 0.18 at a cut, because it saturates, while NDMI falls 0.25–0.40. NDMI
is also far less sensitive to residual haze, because both of its bands sit at
long wavelengths.

### 5.2 The rules, in the order the data forced them

| rule | why it exists |
|---|---|
| drop of at least 0.20, and at least 35 % of the field-year's amplitude | a cut is a real fall, not noise |
| drop **rate** at least 0.010 / day | mowing is a step and drought senescence a slope (a cut drops about 0.25 in 5 days; senescence about 0.30 in 60). Without this rule the detector walks down a decline and calls each stage a cut |
| NDVI must fall in proportion | a real cut removes biomass |
| the next observation within 12 days must stay down | cut grass stays short for weeks, so an instant bounce-back is atmosphere |
| an observation must follow within 20 days | an event with nothing after it cannot be verified, and the haze filter is least reliable at the end of a series |
| a trough with a steep recovery also counts | a cut inside a data gap has no visible falling limb, but its regrowth is still there |
| … only if the canopy was high earlier | spring green-up is also a low point followed by a steep rise; without this guard every field gained a spurious March cut |
| events closer than 21 days merge | the agronomic minimum interval between cuts |

An event is dated at the first observation that is clearly down. That point can
be observed directly, so a human labelling by eye and the detector follow the
same rule. It sits a median of about 3 days after the actual cut, equally in
every year.

The cut series are extracted at a Cloud Score+ threshold of 0.50, looser than
the 0.60 used for condition. The detector was labelled, scored and frozen on
series extracted at 0.50, so changing its input would invalidate the score it
is quoted with. Cut timing also needs every usable date, and the bounce-back
rule above rejects a single hazy one. The cost is that the cut series and the
condition series are not built from exactly the same observations.

### 5.3 Validation, as a supervised-learning problem without a model

- **Reference data.** 65 field-years were labelled by hand in a purpose-built
  widget (`labelling.py`: plain HTML and JS embedded in the notebook, no
  dependencies). The labeller hides the detector's proposals, labels blind, then
  reveals them to compare. A label set copied from the detector measures nothing.
  Labels are kept per site, because a detector validated on Linthebene has not
  been validated on Rheintal.
- **Scoring.** Each detection is matched one-to-one to an unmatched label within
  12 days, so two detections near one cut count as one hit and one false
  positive, not two hits.
- **Train, then freeze.** The rules were tuned on 42 labelled field-years and
  then frozen. Validated on all 65, F1 moved from 0.911 to 0.912, so it did not
  overfit. Re-tuning on the full set would give up that out-of-sample property,
  so the configuration is frozen in `config.yaml`. `evaluate.assert_wired()`
  refuses to run a comparison if a configuration parameter is not actually read
  by the detector. That exact situation once produced a confident and wrong "this
  rule does not help".
- **What the headline needs is stability, not F1.** The claim compares *years*,
  so the threat is an error rate that differs by year. Per-year recall was
  0.812–1.000, with 2026 inside that range.
- **Free ground truth.** Extensive meadows in the valley and hill zones may not
  legally be cut before 15 June, so any earlier detection there is a known
  error. A regulation serves as ground truth.

Precision 0.902, recall 0.922, **F1 0.912**. These scores come from extractions
on the **previous 20 m mask**; the re-run on the 10 m mask is pending
([§5.6](#56-status-and-what-the-full-season-re-run-has-to-watch)).

No per-field cut data ships with this repository, neither real nor anonymised.
The scores above are therefore stated results. Anyone can run the detector and
its scoring harness end to end on synthetic field series (§7), but the real
F1 cannot be reproduced from the repository.

What is *not* measured: truth. The labels are read off the same NDMI curves the
detector uses, so a cut that left no trace in NDMI is invisible to both. This is
agreement with an expert reading of the signal, not a farm diary.

### 5.4 The rejected alternative

A penalised piecewise-linear changepoint fit (`changepoint.py`) was built as the
principled alternative and scored on the same labels. It reached **F1 0.871**,
against the heuristic's 0.911 at the time, mainly by over-segmenting: 32 false
positives against 12. It stays in the code as the record of a negative result.
Two individual rules were also tested and dropped or kept on evidence, not
taste. A minimum-concentration rule *hurt* (F1 0.872 with it, 0.894 without),
because it rejects real cuts whose drop spans two observations. The follow-up
rule stays, because removing it costs three false positives.

There is no learned model here, deliberately. The data are 65 labelled
field-years, and the decision being validated is a comparison between years. A
transparent rule set that was frozen, validated per year, and beat its
principled alternative is the more defensible tool at that sample size. The
obvious next experiment, if wanted, is a small classifier on windowed features
from the same labels, scored through the same harness. Whichever wins, the
comparison is the result.

### 5.5 Which cut metric can carry a between-year claim

- **Cut count** is not comparable between years. Detection sensitivity rises
  steeply with the number of clear observations: about 1.3 cuts are found at
  fewer than 15 images, about 2.9 at 35–45. Drought summers are sunny, and 2026
  had roughly twice the images of 2024. `metrics.equalised` thins every
  field-year to the same number of observations before counting. On the earlier
  run, the equalised counts showed **no** detectable difference between years.
- **Last cut of the season** is the most robust. The detector has to find one
  cut per field-year, not all of them, and the observation bias runs the right
  way: the sunnier year has more images, so a late cut is *easier* to find
  there, not harder. Each field is compared with its own baseline mean, which
  removes site, farm and management style in one step.

### 5.6 Status, and what the full-season re-run has to watch

The cut extractions on disk predate the mask rebuild, and are being re-run now
that the season has closed. The earlier results, all on the previous mask:

- **Last cut about 20 days earlier in 2026** on 12 paired fields, 10 of the 12
  stopping sooner. This comes from hand labels, so it depends on no algorithm;
  the detector corroborates it.
- **No difference in cut count** once observation density is matched.
- **Retracted:** a "share of cuts after DOY 190" finding. Baseline years ran to
  DOY 290 while 2026 stopped at 228, so every September cut counted as "late"
  against a year that could not have one. All years are now truncated at the
  current year's last observation (`detect.truncate_to_current`).

The synthetic data (§7) cannot test how well the detector does on real fields,
but it can expose mechanisms. Using synthetic fields with a planted 2026 shift
and the full 2026 season available (three seeds, two sites):

| comparison ends at DOY | planted last-cut shift | detected (paired median) |
|---|---|---|
| 228 (mid-August) | −33 days | −18 days |
| 248 | −53 days | −31 days |
| 270 | −68 days | +4 days |
| 290 (season end) | −79 days | −11 days |

1. **Up to mid-August the detector follows the planted cuts.** Mean last-cut
   dates per year are within about 4 days of the planted ones. The paired-median
   shift still comes out smaller than planted, because every miss or false
   positive in a single field-year shrinks a pair. So the size of that statistic
   should be quoted cautiously, though its sign holds.
2. **Past mid-August, the regrowth rule fires on the drought recovery.** A deep
   trough followed by a steep rise is exactly the regrowth signature from §5.2.
   On the full season, the frozen rules placed a spurious "cut" at the August
   trough in 41 of about 56 synthetic 2026 field-years. That moved the mean 2026
   last cut from DOY 192 to 232, and the shift collapsed.

The size of the second effect on real fields is unknown. The mechanism,
though, belongs to the frozen rules, not to the synthetic data. So before any
full-season last-cut number is reported:

- fix the comparison end date for all years *before* looking at the result,
  ending before the recovery (mid-August, as the earlier result used);
- count the 2026 detections dated in the trough (about DOY 220–245) and inspect
  them in the labelling widget.

### 5.7 Elevation

Cutting frequency falls with altitude, from 5–6 cuts in the lowlands to 1–2 near
the summering limit. Elevation is carried with every field (SRTM 90 m), so cut
counts are never pooled across altitude unstratified.

---

## 6. Bugs found by checks

None of these was found by luck. Each was caught by a check that should hold
whatever the data, or by a contradiction in the output.

| # | bug | how it was caught | consequence and fix |
|---|---|---|---|
| 1 | Persistence was applied **per cell**, multiplying purity by a "grassland every year" layer | Observed area at purity > 0 came out *smaller* than the polygon area, which is impossible, since the response reaches past every boundary | It zeroed the halo outside fields, and the sample was a shifting set of cells that belonged to no defined set of fields. Rebuilt polygon-first: persistence selects parcels and never touches the raster ([§2.1](#21-which-fields-the-master-parcel-set)). The invariant is now enforced by the demo |
| 2 | Parcels were matched between years by centroid | 56 % of fields appeared to fail persistence, far more land-use change than is plausible | The cause was cadastral redrawing, splitting and merging. Area-overlap matching raised retention from 43 % to 86 % |
| 3 | An index-alignment bug in the overlap routine silently dropped rows whenever its input was a subset | Per-category matching returned 36 % while the all-grassland case passed | The frame built from a dict got a fresh RangeIndex while the geometry kept its own. It is now built from the frame itself |
| 4 | An exception guard made the overlay 100× slower | A three-hour run for an overlay that should take seconds | One malformed geometry sent all ~250,000 pairs to a Python fallback loop. The guard now works per chunk of 5,000 |
| 5 | A figure rebuilt the whole parcel set to compute four sums that were already printed | Eight yearly overlays for one caption | The figure takes the four numbers as arguments |
| 6 | Stale files outlived the method that made them | Slides carried numbers from a superseded mask; a CSV proposed for publication used categories the pipeline no longer emits | The schema was right and the vintage was wrong, and nothing in the code could tell the difference. The fix is structural: provenance records, and a loader that refuses a mismatch ([§8](#8-open-items)) |

Found while packaging this repository:

| # | issue | how it was caught | status |
|---|---|---|---|
| 7 | The mask was built with a single Gaussian (σ 5.30 m) while the deck and `figure_psf` describe box ⊛ Gaussian (σ 4.35 m) | Reading the code that built the mask against the derivation on the slide | Default switched; difference measured (§2.3); published numbers labelled with the model that made them |
| 8 | In the copy of the code handed over for packaging, `trajectory.py` still pooled intensive and extensive meadows | That copy could not produce the published figure's filename: pooling renames the regime before the file is named | One flag; now `analysis.merge_regimes`, off by default. Bug 6 again, caught by the same kind of check |
| 9 | Tiled purity differed from single-pass purity by up to 0.15 | A test that the two must be equal | PIL's polygon fill places edge pixels differently depending on the window offset. Rasterisation moved to GDAL (rasterio), and the results are now identical |
| 10 | "Purity ≥ 1" selected nothing | A threshold panel came out empty | float32 convolution leaves interior cells at 0.99999994. Purity is rounded to six decimals |
| 11 | The 0.25 ha size filter also ran on **past** deliveries | The synthetic cadastre knows which fields are stable, and the persistence test still dropped some | A field whose redrawn outline dipped below 0.25 ha in one year, or which was split into parts below it, counted as "not grassland" that year. On the synthetic cadastre the old rule dropped 520 of 4,040 stable fields (all small) and retained 75.3 %. **Removed**, not patched: the filter duplicated purity's job ([§2.1](#21-which-fields-the-master-parcel-set)). With only the 0.05 ha sliver floor, the funnel recovers exactly what was planted (0 stable fields lost, 85.7 %), and the demo now fails if it does not. `sliver_floor_m2: 2500` reproduces the old behaviour |
| 12 | Geometry repair silently halved figure-eight outlines | A regression test for bug 4's malformed rings, extended to keep area | `buffer(0)` ran first and "repaired" a self-crossing outline into a *valid* polygon with one lobe missing, so `make_valid`, which keeps both, never saw it. Order reversed: `make_valid` first, `buffer(0)` only as a fallback. The synthetic cadastre has no broken outlines, so the demo is unchanged; on the real deliveries the repair count is printed on every build |

---

Each of these has a regression test that reproduces it. Where possible, the
test also runs the old logic and shows that it fails the same check, so the
test is shown to catch the bug, not merely to pass on the fixed code:

| test | bug |
|---|---|
| `test_regressions.py::test_area_at_purity_zero_exceeds_polygon_area` | 1: halo zeroed; the old per-cell multiplication gives a ratio below 1 |
| `test_regressions.py::test_overlap_matching_on_a_subset_keeps_every_row` | 3: the old dict-built frame loses half the rows of a subset |
| `test_regressions.py::test_geometry_repair_survives_pathological_rings` | 4 and 12: self-touching rings and figure-eights, with area kept |
| `test_regressions.py::test_detector_*` | the drop-rate and persistence rules (§5.2): each test fires when its rule is switched off |
| `test_regressions.py::test_config_round_trip` | 6: provenance must describe the configuration that actually ran |
| `test_master_parcels.py` | 11: a field redrawn and split below 0.25 ha keeps its history |
| `test_psf.py::test_implemented_response_has_the_stated_width` | 7: the response the code applies, measured |
| `test_psf.py::test_tiled_purity_equals_single_pass` | 9 |
| `test_cli.py::test_override_after_import_is_refused` | the stale-setting failure behind `assert_wired` (§5.3) |

### How bug 6 became a guarantee

Bug 6 was not fixed by being careful. It was fixed by structure
(`drought.provenance`):

- **Every file carries a record.** It names the settings that shape the
  file, the mask's content hash, resolution, bands and dtype, a hash of the
  file itself, a summary of its content, the package version and a
  timestamp. Earth Engine cannot write such a record, so an export is
  *registered* on arrival.
- **Registration and loading refuse what cannot be current,** and say why.
  That covers a record from other settings (`mask.purity_min differs —
  recorded 0.8, current 0.9`), another mask build, a file changed after it
  was recorded, a missing record, management categories the pipeline does
  not use, missing years, a 20 m or uint8 mask, and a `.geo` column.
- **The Earth Engine scripts take their shared settings from `config.yaml`**
  through `drought gee render`. Each rendered copy writes a hash of its
  settings and mask into its export file names, so registration can tell
  which configuration produced a download.
- **`output/results.json` never mixes vintages.** Each section is stamped with
  the settings and mask it was computed under. A stale section is dropped on
  the next write, and `results.json` refuses to load while one remains.

The CSV that was nearly published, with `mown`/`grazed` categories, a 2026
series ending in mid-August, and a 20 m uint8 mask beside it, would now be
stopped at registration. `tests/test_provenance.py` covers each case.

## 7. The synthetic data, and what it can and cannot show

The real inputs cannot be published, so `drought.synthetic` generates stand-ins
in the same files and schemas, deterministically from one seed:

- **cadastre**: nine annual deliveries with redrawn outlines, splits, merges,
  object IDs that are not stable between years, category changes, and code 930
  in the latest delivery only;
- **trajectories**: the depth and onset of the 2026 collapse are solved so that
  the trajectory module's own `season_integral` returns the real result;
- **field series**: written in Earth Engine export format, with clouds, partial
  coverage, haze, single-date dips, mowing, drought decline and stalled
  regrowth; plus simulated hand labels;
- **station series**: written in MeteoSwiss format.

`make demo` prints each synthetic number beside the real one:

| quantity | synthetic | real | |
|---|---|---|---|
| intensive median ha / < 0.5 ha / shape | 1.11 / 21 % / 5.72 | 1.11 / 21 % / 5.72 | calibrated |
| extensive median ha / < 0.5 ha / shape | 0.45 / 58 % / 6.29 | 0.45 / 58 % / 6.29 | calibrated |
| mown median ha / < 0.5 ha / shape | 0.93 / 26 % / 5.80 | 0.93 / 26 % / 5.80 | emergent |
| integrals (four) | −33 / −35 / −23 / −30 %, rank 1/8 | same | calibrated |
| precipitation, P − ET0 to 10 Aug | 501 mm, −84 mm | same | calibrated |
| area ratio intensive, > 0 / ≥ 0.9 | 1.41 / 0.78 | 1.54 / 0.84 | emergent |
| area ratio extensive, > 0 / ≥ 0.9 | 2.73 / 0.48 | 3.40 / 0.56 | emergent |
| master set retained | 85.7 % | 86.0 % (old filter, provisional) | emergent |
| detector F1, heuristic / changepoint | 0.87–0.90 / 0.80–0.83 | 0.91 / 0.87 | emergent |

The synthetic cadastre contains no fields below 0.25 ha, because the real
geometry targets were measured on fields above that size. So it shows the effect
of removing the filter on *persistence*, but not on the *parcel count*. On the
real data the count will grow by the grassland fields between 0.05 and 0.25 ha,
a number nobody has measured yet. Purity then decides which of them contribute
cells.

*Calibrated* means the generator was solved to reproduce the real value, so a
match only shows the pipeline is wired correctly. *Emergent* means the method
produced it from synthetic input. The area ratios land close to the real ones
without being targeted, which is mild evidence that the synthetic landscape is
realistic where it matters for purity.

**The detector numbers are not evidence about the detector.** On the first
synthetic series, regrowth set in immediately after each cut and the frozen
detector found only half the cuts. The series generator was changed, not the
detector. Cut grass now stays short for about a week and then regrows over
about three weeks, as the detector's own documentation describes real fields.
The detector then scored close to its real F1. That makes the synthetic F1 a
check that the synthetic data look like the data the detector was built for. It
says nothing about the detector's quality on real fields.

What the synthetic data cannot do at all: show that the drought happened. Every
real result above stands or falls on the real data.

---

## 8. Open items

- **Citation for the Sentinel-2 system MTF of 0.25 at Nyquist.** The whole of
  §2.3 rests on it.
- **Purity sensitivity on the 10 m mask** (§2.5). The script is ready: every
  parcel admitted, the Alpstein summering window, Cloud Score+ 0.60.
- **Cut extraction on the 10 m mask, full season** (§5.6), and with it the
  detector re-score.
- **Re-measure the funnel and geometry** on the real data without the 0.25 ha
  filter (§2.1). Every number in §2.1, §2.2 and §2.4 is provisional until then.

