# docs/figures/real/

Figures from the **real** run, copied in by hand. Nothing in the repository
generates them, because the inputs behind them are not published.

Rules for anything placed here:

1. **Aggregates only.** Site-mean trajectories, station series, and
   distributions and totals with no geometry and no parcel IDs
   (`figure_parcel_size`, `area_check`). No map, purity patch or parcel
   outline from the cantonal data: those are derived from the Nutzungsflächen
   and are not publishable. Sentinel-2 true-colour imagery
   with no parcel outlines is open Copernicus data and is allowed.
2. **Current pipeline only.** A figure is accepted only if it was produced with
   `analysis.merge_regimes: false`, and (for the mask figures) without the old
   0.25 ha filter.
3. **Captioned with its run date** in `docs/gallery.md` and the README.

Expected files:

| file | from | run date |
|---|---|---|
| `vi_trajectory_NDMI_intensive_linthebene.png` | `drought trajectory plot --site linthebene` | 14 Oct 2026 |
| `vi_trajectory_NDMI_extensive_linthebene.png` | same | 14 Oct 2026 |
| `vi_trajectory_NDMI_intensive_rheintal.png` | `drought trajectory plot --site rheintal` | 14 Oct 2026 |
| `vi_trajectory_NDMI_extensive_rheintal.png` | same | 14 Oct 2026 |
| `precip_stg_2026.png` | `drought precip --download` | 14 Oct 2026 (data to 12 Oct) |
| `parcel_size.png` | `drought mask figures` | 14 Oct 2026 |
| `area_check_all.png` | same | 14 Oct 2026 |
| `mask_overview_intensive.png` | `drought ee rgb`, then `drought mask overview`. Shows field outlines: the author's explicit exception to rule 1 | 17 Oct 2026 |
| `overview.png` | README figure: upper panels of `precip_stg_2026.png` and `vi_trajectory_NDMI_intensive_linthebene.png` | 20 Oct 2026 |
| `field_anomaly_and_mowing.png` | README figure: per-field NDMI z-score map (slide 22, Code Editor, presentation mask) beside `mowing_one_field_2025_2026.png`. Shows field outlines: the author's explicit exception to rule 1 | Aug 2026 |
| `mixed_pixels_linthebene.png` | README figure: first row of `mask_overview_intensive.png` | 17 Oct 2026 |
| `mowing_one_field_2025_2026.png` | one field's NDMI/NDVI series from the presentation (no ID, no location; the author's explicit exception to rule 1) | Aug 2026 |
| `rheintal_s2_truecolour.jpg` | Earth Engine Code Editor, from the presentation | Aug 2026 (window 15 Jun – 31 Jul) |

`tests/test_docs.py` fails while any image the docs reference is missing, so
the repository cannot be submitted with an empty slot.
