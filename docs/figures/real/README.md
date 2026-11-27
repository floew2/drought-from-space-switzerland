# docs/figures/real/

Figures from the **real** run, copied in by hand. Nothing in the repository
generates them, because the inputs behind them are not published.

Rules for anything placed here:

1. **No parcel data.** Aggregates (site-mean trajectories, station series,
   distributions and totals), Sentinel-2 imagery, and maps that show field
   outlines or per-field values are allowed; the Canton of St. Gallen has
   permitted their use. Parcel IDs and the parcel data itself are never
   published.
2. **Current pipeline only.** A figure is accepted only if it was produced with
   `analysis.merge_regimes: false`, and (for the mask figures) without the old
   0.25 ha filter.
3. **Captioned with its run date** in `docs/gallery.md` and the README.

Expected files:

| file | from | run date |
|---|---|---|
| `vi_trajectory_NDMI_{intensive,extensive}_{linthebene,rheintal,alpstein}.png` | `drought report --only condition` | 27 Nov 2026 |
| `field_anomaly_3sites.jpg` | the three `anomaly_NDMI_2026.png` of `drought report --only anomaly`, side by side | 27 Nov 2026 |
| `vhi_2026.png` | `drought report --only vhi` | 27 Nov 2026 |
| `precip_stg_2026.png` | `drought report --only precip` | 27 Nov 2026 (data to 25 Nov) |
| `overview.png` | README figure: upper panels of `precip_stg_2026.png` and `vi_trajectory_NDMI_intensive_linthebene.png`, with `mowing_one_field_2025_2026.png` below | 27 Nov 2026 |
| `parcel_size.png` | `drought mask figures` | 14 Oct 2026 |
| `area_check_all.png` | same | 14 Oct 2026 |
| `mask_overview_intensive.png` | `drought mask overview` (false colour from the Earth Engine run) | 17 Oct 2026 |
| `mowing_one_field_2025_2026.png` | one field's NDMI/NDVI series from the presentation (no ID, no location) | Aug 2026 |
| `rheintal_s2_truecolour.jpg` | Earth Engine Code Editor, from the presentation | Aug 2026 (window 15 Jun – 31 Jul) |

`tests/test_docs.py` fails while any image the docs reference is missing, so
the repository cannot be submitted with an empty slot.
