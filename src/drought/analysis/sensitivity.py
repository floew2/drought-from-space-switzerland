"""Purity-threshold sensitivity: the anomaly at every threshold, from one pass.

Input: <paths.series>/<site>/sensitivity_<regime>.csv (products/sensitivity.py):
mean index and cell count per purity class per year. The mean at threshold t is
the cell-weighted combination of all classes at or above t — arithmetically
identical to re-extracting at threshold t, but from one pass instead of one per
threshold.
"""

from __future__ import annotations

import pandas as pd

from drought.core.config import CFG, series_dir

__all__ = ["load", "sweep"]


def load(site: str, regime: str) -> pd.DataFrame:
    """A site's sensitivity series for one system, verified."""
    from drought.core.provenance import verify

    f = series_dir(site) / f"sensitivity_{regime}.csv"
    if not f.exists():
        raise FileNotFoundError(f"no {f}; run `drought sensitivity --aoi {site}`")
    verify(f, "sensitivity")
    d = pd.read_csv(f)
    for c in ("year", "pclass", "purity_lo", "mean", "n_px"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    return d.dropna(subset=["mean", "n_px"])


def sweep(d: pd.DataFrame, min_px: int = 200) -> pd.DataFrame:
    """Anomaly of the current year against the baseline median, per threshold.

    Returns:
        threshold, current, baseline, anomaly, n_px; thresholds with fewer
        than `min_px` current-year cells are dropped.
    """
    cur = int(CFG.project.current_year)
    rows = []
    for t in sorted(d.purity_lo.unique()):
        sub = d[d.purity_lo >= t - 1e-9]
        per_year = (sub.assign(w=sub["mean"] * sub.n_px)
                       .groupby("year")[["w", "n_px"]].sum())
        per_year["v"] = per_year.w / per_year.n_px
        if cur not in per_year.index:
            continue
        base = per_year.drop(index=cur)["v"]
        if base.empty:
            continue
        rows.append({"threshold": t, "current": per_year.loc[cur, "v"],
                     "baseline": base.median(),
                     "anomaly": per_year.loc[cur, "v"] - base.median(),
                     "n_px": per_year.loc[cur, "n_px"]})
    t = pd.DataFrame(rows, columns=["threshold", "current", "baseline", "anomaly", "n_px"])
    return t[t.n_px >= min_px].reset_index(drop=True)
