"""Calibration targets for the synthetic data — numbers from the REAL runs.

Every value here was measured on the cantonal data, which cannot be published.
They are the only route by which real data shapes the synthetic data, so each
carries its provenance. Nothing in the pipeline reads this module; only the
generator (to calibrate) and the demo report (to compare) do.

Provenance keys:
    master    `mask.master_parcels()` on the 2017–2025 Nutzungsflächen, full
              extent; box ⊛ Gaussian mask, 0.05 ha floor on 2025 only (2026-10-14)
    geometry  `mask_figures.figure_parcel_size()` on the master set, full extent
    area      `mask.area_table()`, clipped to the raster window
    traj      `trajectory` season NDMI integral, 2026 against the 2019–2025 median
    precip    `precip.main()`, MeteoSwiss station St. Gallen, to 12 October 2026
    detect    `evaluate.score()` on 65 hand-labelled parcel-years — extractions on
              the PREVIOUS 20 m mask, not yet re-run on the 10 m mask
"""

from __future__ import annotations

# --- master parcel set (provenance: master) ----------------------------------
# Measured 2026-10-14: no size filter, the 0.05 ha sliver floor on the 2025
# delivery only (METHOD §2.1, §6 #14). The run with the 0.25 ha filter on every
# delivery (#11) had 40,205 / 34,586 / 80,850 ha.
FUNNEL = {
    "grassland_2025": 66_843,     # grassland parcels >= 0.05 ha in 2025
    "dropped_not_grassland": 7_133,
    "dropped_changed_category": 5_378,
    "exempt_930": 2_418,          # Sömmerung, present only in the latest delivery
    "remain": 54_332,             # 81.3 %
    "remain_ha": 83_662,
}
HELD_ALL_YEARS_SHARE = 0.777      # of the 2025 grassland parcels, 9 of 9 years

CATEGORY_TOTALS = {               # parcels, ha — master set, full extent
    "mown": (45_331, 47_044),
    "grazed": (9_001, 36_619),
    "intensive": (32_991, 43_359),
    "extensive": (12_332, 3_676),
}

# --- parcel geometry (provenance: geometry) -----------------------------------
# median ha, share < 0.5 ha, median shape index (perimeter / sqrt(area))
GEOMETRY = {
    "intensive": {"median_ha": 0.67, "lt_05": 0.43, "median_shape": 5.56},
    "extensive": {"median_ha": 0.19, "lt_05": 0.85, "median_shape": 6.44},
    "mown": {"median_ha": 0.42, "lt_05": 0.54, "median_shape": 5.75},
    "grazed": {"median_ha": 0.58, "lt_05": 0.46, "median_shape": 6.35},
}
# The shape index has a long upper tail: 2,344 of 99,655 parcel-category
# entries (2.4 %) exceed 16. Under the old 0.25 ha filter it was about two in
# 5,000. The small parcels the filter removed include many narrow strips.
SHAPE_TAIL_GT16 = 2_344 / 99_655

# --- area retained by purity threshold (provenance: area) ----------------------
# observed area / polygon area at purity > 0 and at >= 0.9, both clipped to the
# raster window. These are NOT calibration targets: they emerge from geometry +
# PSF, and the demo reports how close the synthetic landscape comes. The ratio
# at > 0 must exceed 1.
AREA_RATIO = {
    "mown": (1.48, 0.80),
    "grazed": (1.45, 0.87),
    "intensive": (1.47, 0.80),
    "extensive": (3.62, 0.45),
}

# --- grassland condition (provenance: traj) -----------------------------------
# season NDMI integral, 2026 vs the 2019–2025 median, %; rank 1 of 8 everywhere.
# Run of 2026-10-14, DOY 121–273 (the preliminary run had -33/-35/-23/-30 to
# DOY 245–250; METHOD §3.4).
INTEGRAL_PCT = {
    ("linthebene", "intensive"): -27.3,
    ("linthebene", "extensive"): -28.9,
    ("rheintal", "intensive"): -22.7,
    ("rheintal", "extensive"): -28.2,
}
NDMI_MIN_AUGUST = (0.09, 0.15)    # the 2026 minimum, mid-August
INTEGRAL_END_DOY = (273, 273)     # real integrals run to here (the window end)
SHARED_DAYS = (103, 114)          # days observed in every year, of 153
# The synthetic 2026 season still ends in early September (END_2026 in
# series.py), as the preliminary exports did, and its draws are accepted at
# 75-100 shared days. The integrals are solved to the values above regardless;
# extending the synthetic season to October is not needed to test the method.
SYNTHETIC_SHARED_DAYS = (75, 100)

# --- precipitation (provenance: precip) ---------------------------------------
PRECIP = {
    "end_doy": 285,               # 12 October 2026
    "p_2026_mm": 690.0,
    "p_median_mm": 1135.0,        # 61 % of median, rank 1 of 31
    "wb_2026_mm": -74.0,          # P - ET0 (Hargreaves)
    "wb_median_mm": 494.0,
}

# --- cut detector (provenance: detect) ----------------------------------------
DETECTOR = {"labelled_parcel_years": 136, "precision": 0.841, "recall": 0.898,
            "f1": 0.868, "f1_without_current": 0.897, "precision_current": 0.642,
            "changepoint_f1": 0.805}
