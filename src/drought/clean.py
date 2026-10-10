"""Cleaning: per-parcel coverage and residual haze."""

from __future__ import annotations

import pandas as pd

from .config import CFG


def clean(df: pd.DataFrame, cfg=None, verbose: bool = True) -> pd.DataFrame:
    """Coverage filter, then the haze filter.

    Haze rows are FLAGGED, not dropped: detection excludes them, but the
    labelling widget and the diagnostic plots show them, so you can see what
    was thrown away and judge whether it deserved to be.
    """
    c = cfg or CFG.clean
    df = df.copy()

    # Coverage is judged PER PARCEL. A 2 ha field does not care whether 80% of
    # the wider site was clear, so the site-level rule used for the trajectory
    # work would discard usable data here.
    best = df.groupby("parcel_id")["n_px"].transform("max")
    if best.max() <= 1:
        if verbose:
            print("  [!] no pixel counts in this export — coverage filter "
                  "skipped")
    else:
        before = len(df)
        df = df[df["n_px"] / best >= c.min_parcel_coverage].copy()
        if verbose:
            print(f"  coverage filter removed {before - len(df)} of {before} "
                  f"parcel-dates below {c.min_parcel_coverage:.0%} of the "
                  "parcel's own best observation")

    # Haze looks exactly like a cut: it depresses NDVI while barely touching
    # NDMI, because scattering inflates short wavelengths and leaves SWIR
    # almost alone. A real cut removes biomass and drops BOTH. Use the
    # disagreement, so no absolute threshold has to be guessed.
    def mark(g):
        g = g.sort_values("doy").copy()
        z = {}
        for ix in ("NDVI", "NDMI"):
            res = g[ix] - g[ix].rolling(c.haze_window, center=True,
                                        min_periods=2).median()
            mad = (res - res.median()).abs().median() * 1.4826
            z[ix] = res / max(mad, 1e-3)
        g["haze"] = (z["NDVI"] < c.haze_ndvi_z) & \
                    ((z["NDMI"] - z["NDVI"]) > c.haze_gap)
        return g

    df = (df.groupby(["parcel_id", "year"], group_keys=False)
            .apply(mark, include_groups=False)
            .join(df[["parcel_id", "year"]]))
    if verbose:
        n = int(df["haze"].sum())
        print(f"  haze filter flagged {n} of {len(df)} parcel-dates "
              f"({100 * n / len(df):.1f}%)")
    return df
