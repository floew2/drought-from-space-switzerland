# docs/figures/real/

Figures from the **real** run, copied in by hand. Nothing in the repository
generates them, because the inputs behind them are not published.

Rules for anything placed here:

1. **Aggregates only.** Site-mean trajectories, station series, and
   distributions and totals with no geometry and no parcel IDs
   (`figure_parcel_size`, `area_check`). No map, purity patch or parcel
   outline from the cantonal data: those are derived from the Nutzungsflächen
   and are not publishable. The purity map and explainer in the gallery come
   from the synthetic demo for that reason.
2. **Current pipeline only.** A figure is accepted only if it was produced with
   `analysis.merge_regimes: false`, and (for the mask figures) without the old
   0.25 ha filter.
3. **Captioned with its run date** in `docs/gallery.md` and the README.

Expected files:

| file | from | run date |
|---|---|---|
| `vi_trajectory_NDMI_intensive_linthebene.png` | `trajectory.plot(site="linthebene")` | *to be filled in* |
| `vi_trajectory_NDMI_extensive_linthebene.png` | same | |
| `vi_trajectory_NDMI_intensive_rheintal.png` | `trajectory.plot(site="rheintal")` | |
| `vi_trajectory_NDMI_extensive_rheintal.png` | same | |
| `precip_stg_2026.png` | `precip.main()` | |
| `parcel_size.png` | `mask.figure_parcel_size(master=...)` | |
| `area_check_all.png` | `mask.area_check(master=...)` | |

`tests/test_docs.py` fails while any image the docs reference is missing, so
the repository cannot be submitted with an empty slot.
