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
  individual parcel failed. Under the old 0.25 ha filter, centroid matching
  kept 43 % and overlap matching 86 %. On every field down to 0.05 ha, overlap
  matching keeps 81 %.
- **Code 930 (Sömmerung) is exempt** from the persistence test, rather than
  failing it. It appears only in the most recent delivery, so there is nothing
  to check it against. That is a schema change, not a change of land use.
- **Persistence is a property of the parcel set and is never multiplied into
  the raster.** [§6](#6-bugs-found-by-checks), bug 1, explains why this matters.
- **There is no size filter.** Whether a field is large enough for Sentinel-2
  is decided by pixel purity (§2.3–2.4), and by nothing else. A field that
  yields no cell above the purity threshold drops out on its own. The only floor
  is a **sliver floor at 0.05 ha**, applied to the 2025 delivery only. It removes
  cadastral offcuts (digitising leftovers and boundary slivers) from the set of
  fields to analyse. It is named for that job, not for observability. It is
  **not** applied to past deliveries. There a polygon of any size can only add
  to the coverage a 2025 field is credited with, so a floor can only fail a
  stable field whose outline was redrawn or split below it that year
  ([§6](#6-bugs-found-by-checks), #14).

An earlier version also applied a 0.25 ha filter, in every year. That was a
second mechanism doing purity's job, and running it on past deliveries failed
stable fields ([§6](#6-bugs-found-by-checks), #11). Size is also a poor proxy
for what the sensor can resolve. Here is the share of isolated rectangular
fields, at random orientation, that yield at least one cell with purity ≥ 0.9
under the box ⊛ Gaussian response:

| shape index | 0.05 ha | 0.10 | 0.15 | 0.20 | 0.25 | 0.30 | 0.40 | 0.50 |
|---|---|---|---|---|---|---|---|---|
| 4.0 (square) | 45 % | 100 % | 100 % | 100 % | 100 % | 100 % | 100 % | 100 % |
| 5.56 (median intensive) | 0 % | 0 % | 6 % | 98 % | 100 % | 100 % | 100 % | 100 % |
| 6.44 (median extensive) | 0 % | 0 % | 0 % | 0 % | 94 % | 100 % | 100 % | 100 % |
| 8.0 | 0 % | 0 % | 0 % | 0 % | 0 % | 0 % | 66 % | 100 % |

(200 random placements and orientations per cell.)

A compact 0.1 ha field is observable; a 0.3 ha strip with shape index 8 is
not. A fixed area cut gets both wrong. Purity gets both right, by definition.

Real funnel, full extent, 14 October 2026:

```
 66,843  grassland parcels >= 0.05 ha in 2025
  7,133  dropped — not grassland in enough years
  5,378  dropped — changed management category
  2,418  exempt (code 930, latest delivery only)
 54,332  remain (81.3 %), 83,662 ha;  77.7 % held their category all nine years
```

Under the old 0.25 ha filter the same cadastre gave 40,205 parcels, of which
34,586 (86.0 %, 80,850 ha) remained. Without the filter the master set holds
57 % more fields (54,332) but only 3.5 % more area: the fields the filter
removed were small. Whether they
contribute anything is then purity's decision (§2.4).

The analysis keeps two category systems, and they do not nest:

| grouping | LNF codes | the question it answers |
|---|---|---|
| intensive / extensive | 613, 650 / 611, 612, 634, 635, 693, 694 | the farming system: 4–6 cuts from early May, against about 2 cuts with none before 15 June |
| mown / grazed | from `config.yaml` | whether a cut signal exists at all |

Intensive and extensive meadows are never pooled. Averaging two cutting regimes
gives a curve that belongs to neither.

### 2.2 Why purity matters: the fields are small and mostly edge

Real master set, full extent:

| category | parcels | median ha | < 0.5 ha | median shape index | 10 m pixels in a median field |
|---|---|---|---|---|---|
| intensive | 32,991 | 0.67 | 43 % | 5.56 | 67 |
| extensive | 12,332 | 0.19 | 85 % | 6.44 | 19 |
| mown | 45,331 | 0.42 | 54 % | 5.75 | 42 |
| grazed | 9,001 | 0.58 | 46 % | 6.35 | 58 |

The shape index is perimeter / √area, where 4.0 is a perfect square. A median
intensive field holds 67 cells of 10 m; a median extensive meadow holds 19. And
a shape index near 6 means a long, narrow field, much of it within a pixel or
two of its boundary. 2.4 % of the parcels have a shape index above 16: strips
along roads, streams and field edges. Under the 0.25 ha filter it was about
two in 5,000, because most strips are small.

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
Gaussian. The preliminary real numbers came from that single-Gaussian mask;
every real number now in this document comes from box ⊛ Gaussian.

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

Real data, master set clipped to the raster window (box ⊛ Gaussian response):

| category | parcels | polygon ha | observed / polygon area at purity > 0 | at ≥ 0.9 |
|---|---|---|---|---|
| mown | 37,483 | 42,052 | 1.48 | 0.80 |
| grazed | 7,488 | 18,772 | 1.45 | 0.87 |
| intensive | 27,400 | 39,039 | 1.47 | 0.80 |
| extensive | 10,076 | 3,004 | 3.62 | 0.45 |

At the threshold of 0.9, intensive meadows keep 80 % of their mapped area. The
small, convoluted extensive meadows keep 45 %, which is the price of measuring
them without their neighbours mixed in. Both are lower than under the old size
filter (84 % and 56 %): the set now includes the small fields, and those hold
the least clean area. That is the filter's work being done by purity instead,
field by field.

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
it can be swept at the cost of one extraction (`drought ee sensitivity`, then
`drought sensitivity plot`; run of 15 October 2026, 10 m mask). The measure is
the mean 2026 NDMI anomaly against 2019–2025 over the site's window:

| site / system | anomaly at 0.00 | at 0.90 (in use) | at 0.95 | spread 0.50–0.95 | 0.90 vs 0.95 |
|---|---|---|---|---|---|
| Linthebene intensive | −0.117 | −0.122 | −0.122 | 3.2 % | 0.7 % |
| Linthebene extensive | −0.096 | −0.114 | −0.118 | 14.5 % | 3.1 % |
| Rheintal intensive | −0.085 | −0.089 | −0.090 | 3.8 % | 0.9 % |
| Rheintal extensive | −0.077 | −0.086 | −0.088 | 8.8 % | 2.9 % |
| Alpstein intensive | −0.098 | −0.102 | −0.103 | 4.3 % | 1.2 % |
| Alpstein extensive | −0.093 | −0.112 | −0.116 | 14.0 % | 3.3 % |

- **The sign never changes.** At every threshold, at every site, for both
  systems, the 2026 anomaly is negative.
- **The size does change, in one direction.** The anomaly grows steadily as the
  threshold tightens, with no clear plateau, so there is no "safe range". The change is
  three to four times larger for extensive meadows, whose fields are small and
  mostly edge (§2.2). That is the expected signature of mixed pixels: cells
  that are partly road, hedge or arable dilute the meadow's signal.
- **So the threshold in use understates, slightly.** At 0.90 the anomaly is
  within 0.7–3.3 % of the strictest value tested. The reported drought signal
  is therefore, if anything, a little too weak, never inflated by contamination.
- **One comparison is not robust: intensive against extensive at Alpstein.**
  Below a threshold of about 0.5, intensive meadows show the larger anomaly;
  above it, extensive ones do. At Linthebene and Rheintal, intensive meadows
  show the larger anomaly at every threshold. No claim here compares the two
  systems at Alpstein.

What the sweep cannot separate is *why* the anomaly drifts. Two causes are
possible. At high thresholds, fewer mixed pixels remain (contamination), but
small fields also drop out of the sample altogether (composition). A second
sweep restricted to parcels that keep a cell of purity ≥ 0.90
(`gee.strict_purity: 0.9`) would split the two. It has not been run, because
it would change the configuration hash of every other export, and the
practical conclusion above does not depend on the answer.

On the previous 20 m mask, the anomaly moved by under 1 % across the whole
range. That sweep had a design flaw. It admitted only parcels with at least
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
  exactly the days that **all** years observed: 103–114 of the 153 in the window,
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
| intensive | −27 % | −23 % |
| extensive | −29 % | −28 % |

Run of 14 October 2026, DOY 121–273, 2019–2025 baseline.

Every one of the four was the lowest of eight years. Spring ran at or above
normal; the index collapsed through July to 0.09–0.15 in mid-August and
recovered sharply from early September. That is a severe midsummer drought that
broke, not a failed year. The two management systems have very different
cutting regimes, so their agreement rules out mowing as the explanation.
Rheintal is consistently milder, as its position near the Rhine and its
shallower groundwater would predict.

**Robustness.** Each variant below changes one choice. 2026 ranks last of eight
in every cell:

| | final | to DOY 245 | haze filter off |
|---|---|---|---|
| Linthebene intensive | −27.3 % | −24.2 % | −27.6 % |
| Linthebene extensive | −28.9 % | −29.1 % | −26.1 % |
| Rheintal intensive | −22.7 % | −23.7 % | −22.9 % |
| Rheintal extensive | −28.2 % | −29.4 % | −28.3 % |

The haze filter removed more dates in 2026 (80) than in any other year (46–57).
But 2026 also had the most dates, so the rate is 19.7 % against 13.4–18.7 %.
Switching the filter off moves three integrals by under half a point and makes
Linthebene extensive less negative by 2.8 points. That is the direction METHOD
§3.1 expects: haze depresses NDVI rather than NDMI, so removing hazy dates
removes noise, not drought.

**What changed from the preliminary run.** A run on the earlier mask gave −33 %
and −35 % at Linthebene, and −23 % and −30 % at Rheintal. Rheintal is unchanged
within two points. Linthebene is smaller by six points in both systems. The
longer window does not explain it: cut at DOY 245, close to where the
preliminary integrals ended, Linthebene intensive is −24 %. Between the two
runs the response model, the parcel set (no size filter), the date coverage
(the preliminary export ended around DOY 245–250) and the Earth Engine implementation
all changed, and the preliminary export cannot be re-read to separate them. The
numbers above are the ones this repository reproduces. The preliminary ones are
recorded here because they were presented before.

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

To 12 October 2026: **690 mm against a median of 1,135 mm**. That is 61 % of
the median and the driest of 31 years in the comparison window. (To 10 August
it was 501 mm against 826, also 61 %: the deficit did not close in autumn.)

Precipitation alone understates a hot drought, so the second panel is the
climatic water balance P − ET0, with ET0 from Hargreaves. This is an index of
atmospheric demand against supply, **not** a soil water balance. It stood at
**−74 mm against a median of +494 mm** on 12 October. Because the quantity
crosses zero, it is reported as a difference, never as a percentage. "−15 % of
normal" would be a serious misreading of a 568 mm shortfall.

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

- **Reference data.** Field-years were labelled by hand in a purpose-built
  widget (`labelling.py`: plain HTML and JS embedded in the notebook, no
  dependencies). The labeller hides the detector's proposals, labels blind, then
  reveals them to compare. A label set copied from the detector measures nothing.
  All labels are from **Linthebene**. 65 existed when the detector was frozen.
  Labelling continued afterwards, and the set now holds 143 field-years with
  391 cuts. Labels are kept per site, because a detector validated on
  Linthebene has not been validated on Rheintal, and it has not been.
- **Scoring.** Each detection is matched one-to-one to an unmatched label within
  12 days, so two detections near one cut count as one hit and one false
  positive, not two hits.
- **Train, then freeze.** The rules were tuned on 42 labelled field-years and
  then frozen. Validated on all 65, F1 moved from 0.911 to 0.912, so it did not
  overfit. The configuration has not changed since, so the 78 field-years
  labelled later were never seen during tuning; they are a genuine hold-out
  (scored here together with the first 65, not separately). Re-tuning on the
  full set would give up that out-of-sample property, so the configuration is
  frozen in `config.yaml`. `evaluate.assert_wired()`
  refuses to run a comparison if a configuration parameter is not actually read
  by the detector. That exact situation once produced a confident and wrong "this
  rule does not help".
- **What the headline needs is stability, not F1.** The claim compares *years*,
  so the threat is an error rate that differs by year. That threat is real:
  recall holds in 2026, precision does not (table below).
- **Free ground truth.** Extensive meadows in the valley and hill zones may not
  legally be cut before 15 June, so any earlier detection there is a known
  error. A regulation serves as ground truth.

Re-scored on the 10 m mask, full season (DOY 60–305), over the 136 of 143
labelled field-years that could be scored (382 cuts), using
`drought cuts score|errors --labelled`:

| | precision | recall | F1 |
|---|---|---|---|
| all years | 0.841 | 0.898 | **0.868** |
| without 2026 | 0.901 | 0.892 | 0.897 |
| 2026 only (61 hits, 34 false positives, 5 misses) | 0.642 | 0.924 | 0.758 |

The drop is a precision drop in one year. §5.6 shows where it comes from. Of
the 39 misses, 17 are the last cut of their season. (The earlier
F1 0.912 was on 65 labels, the previous 20 m mask, and extractions ending in
mid-August 2026.)

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
positives against 12. On the re-scored 136 field-years it reached F1 0.805,
against 0.868. It stays in the code as the record of a negative result.
Two individual rules were also tested and dropped or kept on evidence, not
taste. A minimum-concentration rule *hurt* (F1 0.872 with it, 0.894 without),
because it rejects real cuts whose drop spans two observations. The follow-up
rule stays, because removing it costs three false positives.

There is no learned model here, deliberately. The data are 143 labelled
field-years from one site, and the decision being validated is a comparison between years. A
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

### 5.6 The drought recovery fools the detector, as predicted

The full per-site cut extraction on the 10 m mask exceeds the free Earth Engine
quota, so only the labelled field-years were re-extracted (§5.3). The earlier
results below are on the previous mask:

- **Last cut about 20 days earlier in 2026** on 12 paired fields, 10 of the 12
  stopping sooner. This comes from hand labels, so it depends on no algorithm;
  the detector corroborates it. *Superseded by the label-only result below.*
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

**The real labels confirm it.** The plan was written down before the real
re-score (`drought cuts errors --labelled`). The 34 false positives of 2026, by date:

| DOY | 180–200 | 200–220 | 220–240 | 240–260 | 260–280 |
|---|---|---|---|---|---|
| 2026 false positives | 2 | 5 | **18** | 5 | 4 |

23 of 34 fall between DOY 220 and 260, the trough and the recovery. Without
2026, precision is 0.901; with it, 0.841. So:

- the detector's last cut in 2026 is not used past DOY 228. The fixed
  end date stays at 228 for every year;
- the last-cut shift is computed from the hand labels alone
  (`drought cuts labels`, `metrics.labels_last_cut_paired`), with the detector
  only as corroboration. The labels are read off the same curves, but a human
  labeller sees that a recovery is not a mowing, and the frozen rules do not.

**The last-cut shift, from the hand labels alone** (`drought cuts labels`,
every year ending at DOY 228, 15 October 2026). 22 Linthebene fields are
labelled in 2026 and in at least one earlier year:

- mean last cut DOY 201.5 in the earlier years, DOY 189.7 in 2026;
- median shift **−11.9 days**; 15 of the 22 fields stopped sooner.

This is weaker evidence than the earlier "20 days" suggested. Most fields have
only one baseline year (median 1), so each field's baseline is a single
season, not an average. A sign test on 15 of 22 gives p ≈ 0.13 (two-sided),
not significant on its own. The labelled fields were chosen for labelling, not
sampled at random. The direction agrees with the drought, and the size is
plausible, but on this sample it is an indication, not a finding. Labelling
more baseline years for the same 22 fields would strengthen it most cheaply.

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
| 11 | The 0.25 ha size filter also ran on **past** deliveries | The synthetic cadastre knows which fields are stable, and the persistence test still dropped some | A field whose redrawn outline dipped below 0.25 ha in one year, or which was split into parts below it, counted as "not grassland" that year. On the synthetic cadastre the old rule dropped 520 of 4,040 stable fields (all small) and retained 75.3 %. **Removed**, not patched: the filter duplicated purity's job ([§2.1](#21-which-fields-the-master-parcel-set)). `sliver_floor_m2=2500, floor_past_years=True` reproduces the old behaviour |
| 12 | Geometry repair silently halved figure-eight outlines | A regression test for bug 4's malformed rings, extended to keep area | `buffer(0)` ran first and "repaired" a self-crossing outline into a *valid* polygon with one lobe missing, so `make_valid`, which keeps both, never saw it. Order reversed: `make_valid` first, `buffer(0)` only as a fallback. The synthetic cadastre has no broken outlines, so the demo is unchanged; on the real deliveries the repair count is printed on every build |
| 13 | The build summary divided **clipped** pixels by **unclipped** polygons | The first real build printed area ratios at > 0 below 1 (grazed 0.74, mown 1.32), which bug 1 established is impossible | With `mask.clip_bbox_lv95` set, parcels reach past the grid. Their outside area entered the denominator while their pixels could not exist; summering pastures, mostly outside the box, were hit hardest. Only the printed summary was wrong: the mask is unchanged (bit-identical), and `area_table`, which the figures and `results.json` use, always clipped. The summary now clips the same way |
| 14 | The 0.05 ha sliver floor also ran on **past** deliveries | Recalibrating the synthetic cadastre to the real size distribution, which has many fields of 0.05–0.25 ha: hundreds of planted-stable fields failed persistence (557 of 6,473 on the current synthetic cadastre) | Bug 11 at a smaller scale. A field redrawn or split below 0.05 ha in one year lost that year's coverage. The floor now applies to the 2025 delivery only, and synthetic losses fall to 164. Those are 32 fields whose 2025 outline came out below the floor, and 132 narrow strips (median shape index 14). Overlap matching cannot follow a strip a few metres wide through a 0.5 m redraw. **None of the 164 can hold a clean cell on its own** (highest own purity 0.74), so the analysed signal is unaffected; they leave only halo. The demo now fails if a lost stable field could hold a clean cell. On the real cadastre the fix kept 252 more fields (54,080 → 54,332; 80.9 % → 81.3 %), far fewer than the synthetic estimate suggested, so the synthetic redraw model (0.5 m per vertex, every year) seems to move outlines more than real deliveries do |

---

Each of these has a regression test that reproduces it. Where possible, the
test also runs the old logic and shows that it fails the same check, so the
test is shown to catch the bug, not merely to pass on the fixed code:

| test | bug |
|---|---|
| `test_regressions.py::test_area_at_purity_zero_exceeds_polygon_area` | 1: halo zeroed; the old per-cell multiplication gives a ratio below 1 |
| `test_regressions.py::test_area_ratio_reference_is_clipped_to_the_raster` | 13: the unclipped denominator gives a ratio below 1 |
| `test_regressions.py::test_overlap_matching_on_a_subset_keeps_every_row` | 3: the old dict-built frame loses half the rows of a subset |
| `test_regressions.py::test_geometry_repair_survives_pathological_rings` | 4 and 12: self-touching rings and figure-eights, with area kept |
| `test_regressions.py::test_detector_*` | the drop-rate and persistence rules (§5.2): each test fires when its rule is switched off |
| `test_regressions.py::test_config_round_trip` | 6: provenance must describe the configuration that actually ran |
| `test_master_parcels.py` | 11: a field redrawn and split below 0.25 ha keeps its history; 14: a split part below 0.05 ha in a past year does not fail the field |
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
- **The Earth Engine exports run from Python** (`drought ee …`), reading
  `config.yaml` directly; there is no hand-edited copy of a setting to drift.
  The mask asset is compared with the local mask cell by cell before it is
  stamped with the mask's content hash, and an export refuses to start against
  an asset with any other stamp. Exports are named with a hash of every
  Earth Engine setting plus the mask, and each download is registered as it
  arrives. The two remaining Code Editor scripts, a map and a hotspot view,
  get their shared settings from `drought gee render`, checked by CI.
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
| intensive median ha / < 0.5 ha / shape | 0.70 / 42 % / 5.54 | 0.67 / 43 % / 5.56 | calibrated |
| extensive median ha / < 0.5 ha / shape | 0.20 / 84 % / 6.35 | 0.19 / 85 % / 6.44 | calibrated |
| mown median ha / < 0.5 ha / shape | 0.45 / 53 % / 5.76 | 0.42 / 54 % / 5.75 | emergent |
| integrals (four) | −27 / −29 / −23 / −28 %, rank 1/8 | same | calibrated |
| precipitation, P − ET0 to 12 Oct | 689 mm, −75 mm | 690 mm, −74 mm | calibrated |
| area ratio intensive, > 0 / ≥ 0.9 | 1.54 / 0.74 | 1.47 / 0.80 | emergent |
| area ratio extensive, > 0 / ≥ 0.9 | 3.18 / 0.36 | 3.62 / 0.45 | emergent |
| master set retained | 79.1 % | 81.3 % | planted |
| detector F1, heuristic / changepoint | 0.87–0.90 / 0.80–0.83 | 0.91 / 0.87 | emergent |

The synthetic set retains less than the real one (79.1 against 81.3 %)
because the generator redraws every outline by 0.5 m in every year, which loses
the strips described below. Real deliveries apparently move outlines less
(bug 14 cost 252 real fields, not the several hundred per 6,500 it cost
synthetically).

The synthetic cadastre is calibrated to the real size distribution down to the
0.05 ha floor, so it has the many small fields the old filter removed. That is
how bug 14 surfaced: a floor on past deliveries fails small fields that were
redrawn or split. Under the corrected rule the funnel recovers what was planted,
with one class of exception that the demo reports and does not fail. Strips a
few metres wide cannot be followed through a redraw by 80 % overlap. None of
them can hold a clean cell, so they cost halo, not signal.

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
- **Cut extraction on the 10 m mask, full season** (§5.6). The labelled
  field-years are done and the detector re-scored; the full per-site sample
  exceeds the free Earth Engine quota.
- **More baseline labels for the 22 paired Linthebene fields** (§5.6): the
  cheapest way to firm up the last-cut shift.
