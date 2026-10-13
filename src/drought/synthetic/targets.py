"""Calibration targets for the synthetic data — numbers from the REAL runs.

Every value here was measured on the cantonal data, which cannot be published.
They are the only route by which real data shapes the synthetic data, so each
carries its provenance. Nothing in the pipeline reads this module; only the
generator (to calibrate) and the demo report (to compare) do.

Provenance keys:
    master    `mask.master_parcels()` on the 2017–2025 Nutzungsflächen, full
              extent; box ⊛ Gaussian mask, 0.05 ha sliver floor (2026-10-13)
    geometry  `mask_figures.figure_parcel_size()` on the master set, full extent
    area      `mask.area_table()`, clipped to the raster window
    traj      `trajectory` season NDMI integral, 2026 against the 2019–2025 median
    precip    `precip.main()`, MeteoSwiss station St. Gallen, to 10 August 2026
    detect    `evaluate.score()` on 65 hand-labelled parcel-years — extractions on
              the PREVIOUS 20 m mask, not yet re-run on the 10 m mask
"""

from __future__ import annotations

# --- master parcel set (provenance: master) ----------------------------------
# First 10 m build, 2026-10-13: no size filter, but the 0.05 ha sliver floor
# still ran on past deliveries (METHOD §6, #14). The 2025 count is unaffected;
# everything measured on the master set (remain, category totals, geometry,
# area ratios) is re-measured on the corrected rule. The run before, with the 0.25 ha filter on every delivery
# (#11), had 40,205 / 34,586 / 80,850 ha.
FUNNEL = {
    "grassland_2025": 66_843,     # grassland parcels >= 0.05 ha in 2025
    "dropped_not_grassland": 7_540,
    "dropped_changed_category": 5_223,
    "exempt_930": 2_418,          # Sömmerung, present only in the latest delivery
    "remain": 54_080,             # 80.9 %
    "remain_ha": 83_633,
}
HELD_ALL_YEARS_SHARE = 0.773      # of the 2025 grassland parcels, 9 of 9 years

CATEGORY_TOTALS = {               # parcels, ha — master set, full extent
    "mown": (45_108, 47_019),
    "grazed": (8_972, 36_614),
    "intensive": (32_899, 43_349),
    "extensive": (12_202, 3_661),
}

# --- parcel geometry (provenance: geometry) -----------------------------------
# median ha, share < 0.5 ha, median shape index (perimeter / sqrt(area))
GEOMETRY = {
    "intensive": {"median_ha": 0.67, "lt_05": 0.42, "median_shape": 5.56},
    "extensive": {"median_ha": 0.19, "lt_05": 0.85, "median_shape": 6.43},
    "mown": {"median_ha": 0.43, "lt_05": 0.54, "median_shape": 5.75},
    "grazed": {"median_ha": 0.58, "lt_05": 0.45, "median_shape": 6.35},
}
# The shape index has a long upper tail: 2,309 of 99,181 parcel-category
# entries (2.3 %) exceed 16. Under the old 0.25 ha filter it was about two in
# 5,000. The small parcels the filter removed include many narrow strips.
SHAPE_TAIL_GT16 = 2_309 / 99_181

# --- area retained by purity threshold (provenance: area) ----------------------
# observed area / polygon area at purity > 0 and at >= 0.9, both clipped to the
# raster window. These are NOT calibration targets: they emerge from geometry +
# PSF, and the demo reports how close the synthetic landscape comes. The ratio
# at > 0 must exceed 1.
AREA_RATIO = {
    "mown": (1.48, 0.80),
    "grazed": (1.45, 0.87),
    "intensive": (1.47, 0.80),
    "extensive": (3.61, 0.45),
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
