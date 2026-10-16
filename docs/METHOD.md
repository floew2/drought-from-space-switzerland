# Method

How the analysis works, and why it is built this way. For results at a glance
and how to run the code, see the [README](../README.md).

1. [Question and data](#1-question-and-data)
2. [The grassland mask](#2-the-grassland-mask)
3. [Grassland condition](#3-grassland-condition)
4. [Precipitation and water balance](#4-precipitation-and-water-balance)
5. [Cut detection](#5-cut-detection)
6. [Bugs found by checks](#6-bugs-found-by-checks)
7. [Synthetic data](#7-synthetic-data)
8. [Limitations and open items](#8-limitations-and-open-items)

---

## 1. Question and data

How hard did the 2026 drought hit grassland in the Canton of St. Gallen, and
did farmers stop mowing earlier? Two valley-floor grassland areas are analysed,
Linthebene and the St. Galler Rheintal, plus one alpine summering area
(Alpstein) for context.

| input | what | used for |
|---|---|---|
| Sentinel-2 L2A, 10 m, 2019–2026 | spectral indices, Cloud Score+ cloud mask | condition and cuts |
| Nutzungsflächen 2017–2025 | yearly field outlines with land-use codes | which pixels are grassland, of which kind |
| MeteoSwiss station St. Gallen, 1991–2026 | daily precipitation and temperature | independent check |

2026 is compared with 2019–2025. That baseline includes the 2022 drought, so
the 2026 anomaly is, if anything, understated. Everything is in CH1903+ / LV95
(EPSG:2056).

![St. Galler Rheintal near Altstätten: Sentinel-2 true colour, 2019–2025 median against 2026](figures/real/rheintal_s2_truecolour.jpg)

*Sentinel-2 true colour, 15 June – 31 July, identical stretch: the median of
2019–2025 (left) against 2026 (right). Meadows browned field by field.*

---

## 2. The grassland mask

Most Swiss meadows are small, so most 10 m pixels over them also see a road, a
hedge or a neighbouring crop. The mask answers one question per pixel: *how much
of what the sensor records here comes from inside one kind of grassland?*

### 2.1 Which fields: the master parcel set

A field is analysed only if it was grassland **of the same management type**
(intensive or extensive) in every year from 2017 to 2025. Outlines are redrawn,
split and merged between yearly deliveries, so fields are matched by **area
overlap** (≥ 80 %), not by ID or centroid. Centroid matching kept 43 % of
fields; overlap matching keeps 81 %.

```
 66,843  grassland parcels ≥ 0.05 ha in 2025
  7,133  dropped: not grassland in every year
  5,378  dropped: changed management type
  2,418  exempt: summering code 930, present only in 2025
 54,332  remain (81.3 %), 83,662 ha
```

There is no size filter. Whether a field is large enough is decided by pixel
purity (§2.3), not by area. A compact 0.1 ha field can hold a clean pixel; a
0.3 ha strip cannot. Intensive and extensive meadows are never pooled: their
cutting regimes differ (4–6 cuts from early May against about 2, none before
15 June), and an average of the two belongs to neither.

### 2.2 Why purity matters: the fields are small and mostly edge

![Parcel size and shape index by management category](figures/real/parcel_size.png)

| category | parcels | median ha | < 0.5 ha | median shape index | 10 m pixels in a median field |
|---|---|---|---|---|---|
| intensive | 32,991 | 0.67 | 43 % | 5.56 | 67 |
| extensive | 12,332 | 0.19 | 85 % | 6.44 | 19 |

The shape index is perimeter / √area; 4 is a square. A median extensive meadow
holds 19 pixels, and with a shape index above 6 most of them lie within a pixel
or two of the edge. This is why the mask uses the native 10 m grid, not 20 m.

### 2.3 The sensor response

A pixel does not see a sharp 10 × 10 m square: it also draws signal from
around it. Purity is therefore defined on the sensor's **response**, not on
area. Parcels are rasterised at 2.5 m and convolved with a model of the
response, and each 10 m cell gets the share of its signal that comes from inside
the right fields.

![Box convolved with Gaussian optics: FWHM 12.7 m](figures/model/psf.png)

The model is a 10 m detector box convolved with Gaussian optics, derived from
the system MTF (modulation transfer function) at the Nyquist frequency:

| quantity | value | how |
|---|---|---|
| system MTF at Nyquist | 0.25 | chosen inside ESA's specified range, 0.15–0.30 for the 10 m bands ([ESA](https://sentinel.esa.int/web/sentinel/technical-guides/sentinel-2-msi/mission-performance)) |
| detector MTF | 0.637 | box one pixel wide: sinc(½) |
| optics MTF | 0.393 | 0.25 / 0.637 |
| σ of the Gaussian | 4.35 m | from MTF(f) = exp(−2π²σ²f²) |
| width of the response | 12.7 m (FWHM) | box ⊛ Gaussian, measured numerically |

A test measures the width the code actually applies, so the implementation
cannot drift from this derivation. The model is nadir-only and isotropic, and
real optics are not exactly Gaussian.

### 2.4 What the threshold costs

Only cells with purity ≥ 0.90 are used. That costs area, most of all for the
small extensive meadows:

![Three landscapes: Sentinel-2 imagery, field outlines, and purity at three thresholds](figures/real/mask_overview_intensive.png)

*Intensive grassland at Linthebene, in the Rheintal and at Alpstein, 500 m
across. From left: Sentinel-2 false colour, the field outlines, and purity at
three thresholds (white pure, grey mixed, black discarded). The percentages are
retained area against field area: above 100 % at a low threshold, the response
reaches past the fields. From the August 2026 presentation, built with the
earlier single-Gaussian response (σ 5.3 m); `drought mask overview` redraws it
from the current mask.*

![Observed area against polygon area, by purity threshold](figures/real/area_check_all.png)

| category | observed / polygon area at purity > 0 | at ≥ 0.90 |
|---|---|---|
| intensive | 1.47 | 0.80 |
| extensive | 3.62 | 0.45 |

The ratio at purity > 0 **must exceed 1**, because the response reaches past
every field boundary. This invariant caught the worst bug of the project
([§6](#6-bugs-found-by-checks), bug 1), and the demo enforces it.

### 2.5 Is the threshold a lever on the answer?

The threshold is applied in Earth Engine, so it can be swept. The 2026 NDMI
anomaly against 2019–2025, by threshold:

| | at 0.00 | at 0.90 (used) | at 0.95 | spread 0.50–0.95 |
|---|---|---|---|---|
| Linthebene intensive | −0.117 | −0.122 | −0.122 | 3.2 % |
| Linthebene extensive | −0.096 | −0.114 | −0.118 | 14.5 % |
| Rheintal intensive | −0.085 | −0.089 | −0.090 | 3.8 % |
| Rheintal extensive | −0.077 | −0.086 | −0.088 | 8.8 % |

The anomaly is negative at every threshold, and grows as the threshold
tightens: mixed pixels dilute the drought signal, most of all on small
extensive fields. At 0.90 it is within 0.7–3.3 % of the strictest value, so the
reported signal is slightly understated, never inflated. (At Alpstein, which of
the two meadow types shows the larger anomaly flips near a threshold of 0.5, so
the two are not compared there.)

---

## 3. Grassland condition

### 3.1 Site-mean trajectories

For each site, management type and acquisition date, Earth Engine returns the
mean index over pixels with purity ≥ 0.90. Before any curve is drawn:

- **Clouds.** A pixel is used only if its Cloud Score+ value is at least 0.60.
  A date is dropped if less than 80 % of the site is clear, because a partly
  clouded scene samples *different fields*, not a random subset.
- **Overlapping tiles.** Two Sentinel-2 tiles seeing the same ground on the
  same day are merged, weighted by pixel count.
- **Residual haze.** Thin haze passes the cloud mask, but depresses red-band
  NDVI far more than SWIR-based NDMI. Dates where NDVI drops while NDMI does not
  are flagged with robust z-scores (about 17 % of dates, evenly across years).

### 3.2 The season integral

The headline is the area under the NDMI curve from 1 May to 30 September
(DOY 121–273), relative to the 2019–2025 median. Every year is integrated over
exactly the days **all** years observed (103–114 of 153), so a cloudier year is
not penalised for missing data. Gaps longer than 30 days are left empty rather
than interpolated.

### 3.3 Why NDMI

NDMI (near-infrared against shortwave infrared) responds to canopy water and,
unlike NDVI, does not saturate over dense meadow. It has no absolute scale, so
it is only read against each site's own history.

### 3.4 Result

![Linthebene, intensive meadows: 2026 NDMI against the 2019–2025 envelope](figures/real/vi_trajectory_NDMI_intensive_linthebene.png)

| season NDMI integral | Linthebene | St. Galler Rheintal |
|---|---|---|
| intensive meadows | −27 % | −23 % |
| extensive meadows | −29 % | −28 % |

2026 was the lowest of eight years in all four series. Spring ran at or above
normal; the index collapsed through July to 0.09–0.15 in mid-August and
recovered sharply from early September. Both management types agree, so mowing
is not the explanation. The Rheintal, near the Rhine with shallower
groundwater, is consistently milder.

**Robustness.** With the haze filter off, three integrals move by less than
half a point and Linthebene extensive by 2.8 points; 2026 still ranks last
everywhere. Ending the window at DOY 245 instead of 273 changes the integrals
by about 3 points or less.

---

## 4. Precipitation and water balance

![Cumulative precipitation and P − ET0, 2026 against 1991–2020](figures/real/precip_stg_2026.png)

The current year and the 1991–2020 reference come from the same station
series, aligned by calendar day, so a difference cannot come from a change of
instrument. To 12 October 2026: **690 mm against a median of 1,135 mm**, 61 %
and the driest of 31 years.

Rain alone understates a hot drought, so the second panel shows the climatic
water balance P − ET0 (Hargreaves): **−74 mm against a median of +494 mm**.
This is atmospheric supply against demand, not soil moisture.

---

## 5. Cut detection

### 5.1 Why per field

A site mean averages hundreds of fields cut on different days, which smooths the
mowing sawtooth away. So cuts are detected per field, on a sample of fields of
at least 0.5 ha, each shrunk inward by 15 m.

![One meadow, 2025 against 2026: NDMI and NDVI with detected cuts](figures/real/mowing_one_field_2025_2026.png)

### 5.2 The rules

A cut is a sudden drop in NDMI, confirmed by NDVI:

| rule | why |
|---|---|
| drop ≥ 0.20 and ≥ 35 % of the season's range | a real fall, not noise |
| drop rate ≥ 0.010 per day | mowing is a step; drought decline is a slope |
| NDVI falls too | a real cut removes biomass |
| still low at the next observation within 12 days | cut grass stays short; an instant rebound is atmosphere |
| a later observation within 20 days | an event at the very end cannot be verified |
| a trough followed by steep regrowth also counts | catches a cut hidden in a data gap |
| events closer than 21 days merge | the agronomic minimum between cuts |

The cut series use a Cloud Score+ threshold of 0.50, the value the detector was
tuned and scored with.

### 5.3 Validation

The rules were tuned on 42 hand-labelled field-years and then **frozen**. Labels
were made in a purpose-built widget that hides the detector's proposals while
labelling. The label set has since grown to 143 field-years at Linthebene
without any change to the rules, so the later labels are a true hold-out. Each
detection is matched one-to-one to a label within 12 days.

| | precision | recall | F1 |
|---|---|---|---|
| all 136 scorable field-years | 0.841 | 0.898 | **0.868** |
| without 2026 | 0.901 | 0.892 | 0.897 |
| 2026 only | 0.642 | 0.924 | 0.758 |

The labels are read off the same curves the detector uses, so this measures
agreement with an expert reading of the signal, not with a farm diary.

### 5.4 The rejected alternative

A penalised changepoint fit was built as the principled alternative and scored
on the same labels. It over-segments and reaches F1 0.805, so it was rejected.
It stays in the code (`changepoint.py`). With 143 labelled field-years from one
site, a transparent, frozen rule set is the more defensible tool than a learned
model.

### 5.5 Which cut metric can compare years

Cut **counts** are not comparable: sunny drought summers give more clear images,
and more images reveal more cuts. The **last cut of the season** is the robust
metric. Each field is compared with its own earlier years, which removes site
and farm effects.

### 5.6 The drought recovery fools the detector

On synthetic data the frozen rules had shown a weak spot: a deep drought trough
followed by a steep recovery looks exactly like a cut followed by regrowth. The
real labels confirm it. Of the detector's 34 false cuts in 2026, 23 fall between
DOY 220 and 260, the trough and the recovery.

| DOY | 180–200 | 200–220 | 220–240 | 240–260 | 260–280 |
|---|---|---|---|---|---|
| 2026 false positives | 2 | 5 | **18** | 5 | 4 |

So the last-cut comparison ends at DOY 228 (mid-August) for every year, and it
is computed **from the hand labels alone**. 22 fields are labelled both in 2026
and in at least one earlier year:

- the last cut came a median **11.9 days earlier** in 2026;
- 15 of the 22 fields stopped sooner.

Most fields have only one earlier year to compare with, and a sign test gives
p ≈ 0.13. The direction agrees with the drought, but on this sample it is an
indication, not a significant result.

### 5.7 Elevation

Cutting frequency falls with altitude, so every field carries its elevation
(SRTM 90 m) and cut counts are never pooled across altitude.

---

## 6. Bugs found by checks

Each bug below was caught by a check that must hold whatever the data, or by a
contradiction in the output, and the code-level ones have a regression test in
`tests/` that the old logic fails.

| # | bug | caught by |
|---|---|---|
| 1 | persistence multiplied into the raster, zeroing the halo around fields | observed area at purity > 0 below the polygon area, which is impossible |
| 2 | fields matched between years by centroid | 56 % apparent land-use change, far too much |
| 3 | overlap matching dropped rows on a subset | per-category results contradicting the all-grassland case |
| 4 | one bad geometry made the overlay 100× slower | a three-hour run |
| 5 | a figure rebuilt the whole parcel set for four sums | eight overlays for one caption |
| 6 | results outlived the settings that made them | numbers from a superseded mask; now refused by provenance records |
| 7 | mask built with a different response model than documented | comparing the code with the derivation; now a test measures the applied width |
| 8 | intensive and extensive meadows pooled in one code copy | output file names that could not match |
| 9 | tiled purity differed from single-pass purity | a test that they must be equal |
| 10 | "purity ≥ 1" selected nothing | an empty panel (float32 rounding) |
| 11 | a 0.25 ha size filter also failed stable fields in past years | synthetic fields known to be stable being dropped |
| 12 | geometry repair silently halved figure-eight outlines | a test that repair keeps area |
| 13 | a summary divided clipped pixels by unclipped polygons | area ratios below 1 again |
| 14 | the 0.05 ha floor also failed stable fields in past years | the same synthetic check as #11 |

**How bug 6 became a guarantee.** Every data file carries a provenance record:
the settings and mask that produced it, and hashes of both. Loading refuses a
file whose record does not match the current configuration, and names the
setting that changed. Earth Engine exports are named with a hash of their
settings and the mask, and the uploaded mask asset is compared with the local
file cell by cell before any export runs. `results.json` refuses to load while
any section is stale.

---

## 7. Synthetic data

The real parcel data cannot be published, so `drought.synthetic` generates
stand-ins with the same files and schemas: yearly cadastres with redrawn, split
and merged outlines; Earth Engine exports with clouds, haze, mowing and drought;
hand labels; and a MeteoSwiss-format station series. The generator is
calibrated so that the same code reproduces the real headline numbers. Matching
a calibrated number shows only that the pipeline is wired correctly. The area
ratios, which were not targeted, land close to the real ones.

What synthetic data cannot do: show that the drought happened. Every result
above stands on the real data.

---

## 8. Limitations and open items

- **The MTF value.** ESA specifies 0.15–0.30; 0.25 is assumed. Rebuilding the
  synthetic mask at 0.15 and 0.30 (`--set mask.psf_mtf_nyquist=…`) would show
  what the choice costs.
- **Mowing sample.** The last-cut result rests on 22 fields at one site, most
  with one baseline year. More baseline labels for the same fields would
  strengthen it most cheaply. The detector has not been validated outside
  Linthebene.
- **Full cut extraction.** A per-site sample of 250 fields over eight years
  exceeds the free Earth Engine quota; only the labelled field-years were
  extracted.
- **NDMI is a proxy.** It tracks canopy water and density, not yield. The
  integral is a relative measure against each site's own history.
