"""Calibration targets for the synthetic data — numbers from the REAL runs.

Every value here was measured on the cantonal data, which cannot be published.
They are the only route by which real data shapes the synthetic data, so each
carries its provenance. Nothing in the pipeline reads this module; only the
generator (to calibrate) and the demo report (to compare) do.

Provenance keys:
    master    `mask.master_parcels()` on the 2017–2025 Nutzungsflächen, full extent
    geometry  `mask_figures.figure_parcel_size()` on the master set, full extent
    area      `mask.area_table()`, clipped to the raster window
    traj      `trajectory` season NDMI integral, 2026 against the 2019–2025 median
    precip    `precip.main()`, MeteoSwiss station St. Gallen, to 10 August 2026
    detect    `evaluate.score()` on 65 hand-labelled parcel-years — extractions on
              the PREVIOUS 20 m mask, not yet re-run on the 10 m mask
"""

from __future__ import annotations

# --- master parcel set (provenance: master) ----------------------------------
# PROVISIONAL: measured with the old 0.25 ha size filter, which also ran on
# past deliveries (METHOD §6, #11). The filter has been replaced by a 0.05 ha
# sliver floor; the real funnel is re-measured on the next real run. The
# generator plants statuses at these shares, so under the new rule the
# synthetic funnel should reproduce them exactly.
FUNNEL = {
    "grassland_2025": 40_205,     # grassland parcels >= 0.25 ha in 2025
    "dropped_not_grassland": 3_592,
    "dropped_changed_category": 2_027,
    "exempt_930": 2_137,          # Sömmerung, present only in the latest delivery
    "remain": 34_586,             # 86.0 %
    "remain_ha": 80_850,
}

CATEGORY_TOTALS = {               # parcels, ha — master set
    "mown": (28_112, 44_628),
    "grazed": (6_474, 36_223),
    "intensive": (23_758, 42_034),
    "extensive": (4_348, 2_585),
}

# --- parcel geometry (provenance: geometry) -----------------------------------
# median ha, share < 0.5 ha, median shape index (perimeter / sqrt(area))
GEOMETRY = {
    "intensive": {"median_ha": 1.11, "lt_05": 0.21, "median_shape": 5.72},
    "extensive": {"median_ha": 0.45, "lt_05": 0.58, "median_shape": 6.29},
    "mown": {"median_ha": 0.93, "lt_05": 0.26, "median_shape": 5.80},
    "grazed": {"median_ha": 0.99, "lt_05": 0.25, "median_shape": 6.61},
}
# The shape index has a long upper tail: a minority of parcels exceed 10, and
# about two in 5,000 exceed 16.
SHAPE_TAIL_GT16 = 2 / 5_000

# --- area retained by purity threshold (provenance: area) ----------------------
# observed area / polygon area at purity > 0 and at >= 0.9. These are NOT
# calibration targets: they emerge from geometry + PSF, and the demo reports
# how close the synthetic landscape comes. The ratio at > 0 must exceed 1.
AREA_RATIO = {
    "mown": (1.54, 0.84),
    "grazed": (1.49, 0.88),
    "intensive": (1.54, 0.84),
    "extensive": (3.40, 0.56),
}

# --- grassland condition (provenance: traj) -----------------------------------
# season NDMI integral, 2026 vs the 2019–2025 median, %; rank 1 of 8 everywhere
INTEGRAL_PCT = {
    ("linthebene", "intensive"): -33.0,
    ("linthebene", "extensive"): -35.0,
    ("rheintal", "intensive"): -23.0,
    ("rheintal", "extensive"): -30.0,
}
NDMI_MIN_AUGUST = (0.09, 0.15)    # the 2026 minimum, mid-August
INTEGRAL_END_DOY = (245, 250)     # integrals run to here
SHARED_DAYS = (75, 100)           # days observed in every year, of 153

# --- precipitation (provenance: precip) ---------------------------------------
PRECIP = {
    "end_doy": 222,               # 10 August 2026
    "p_2026_mm": 501.0,
    "p_median_mm": 826.0,         # 61 % of median, rank 1 of 31
    "wb_2026_mm": -84.0,          # P - ET0 (Hargreaves)
    "wb_median_mm": 324.0,
}

# --- cut detector (provenance: detect) ----------------------------------------
DETECTOR = {"labelled_parcel_years": 65, "precision": 0.902, "recall": 0.922,
            "f1": 0.912, "recall_per_year": (0.812, 1.000),
            "changepoint_f1": 0.871}
