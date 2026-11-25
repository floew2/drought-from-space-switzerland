"""The per-field anomaly map: each parcel coloured by its z-score.

Red is drier than the parcel's own baseline years, blue wetter, in standard
deviations of the baseline; grey parcels have no z-score (too few clear
observations or baseline years).
"""

from __future__ import annotations

import logging

import matplotlib.pyplot as plt
import pandas as pd

from drought.core.config import CFG, figures_dir, series_dir, window

__all__ = ["figure", "load"]

log = logging.getLogger(__name__)
VMAX = 3.0


def load(site: str) -> pd.DataFrame:
    """A site's per-field anomaly table, verified."""
    from drought.core.provenance import verify

    f = series_dir(site) / "anomaly.csv"
    if not f.exists():
        raise FileNotFoundError(f"no {f}; run `drought anomaly --aoi {site}`")
    verify(f, "anomaly")
    return pd.read_csv(f, dtype={"parcel_id": str})


def figure(site: str) -> dict:
    """Map of one site's per-field z-scores; returns summary numbers."""
    from drought.core import aoi
    from drought.products.common import parcels

    t = load(site)
    g = parcels(site).merge(t[["parcel_id", "regime", "z"]], on="parcel_id")
    cur = int(CFG.project.current_year)
    lo, hi = window(site)

    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(7.2, 7.2), dpi=180, constrained_layout=True)
    a = aoi.get(site)
    x0, y0, x1, y1 = a.bounds
    ax.set_facecolor("#0d1117")
    g[g.z.isna()].plot(ax=ax, color="#484f58", linewidth=0)
    shown = g[g.z.notna()]
    if not shown.empty:
        shown.plot(ax=ax, column="z", cmap="RdBu", vmin=-VMAX, vmax=VMAX, linewidth=0)
        sm = plt.cm.ScalarMappable(cmap="RdBu", norm=plt.Normalize(-VMAX, VMAX))
        cb = fig.colorbar(sm, ax=ax, shrink=0.7, pad=0.01)
        cb.set_label(f"{CFG.anomaly.index} {cur} vs. own baseline (SD)", color="white",
                     fontsize=9)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_color("#30363d")
    ax.set_title(f"{a.label}: per-field {CFG.anomaly.index} anomaly, {cur}\n"
                 f"DOY {lo}-{hi}, against {CFG.project.baseline_years[0]}-"
                 f"{CFG.project.baseline_years[-1]}", fontsize=11, fontweight="bold",
                 color="white", loc="left")
    ax.set_xlabel("Sentinel-2, swissEO S2-SR (swisstopo); parcels: Nutzungsflächen, "
                  "Kanton St.Gallen. Grey: no z-score.", fontsize=7, color="#6e7681",
                  loc="left")
    out = figures_dir(site) / f"anomaly_{CFG.anomaly.index}_{cur}.png"
    fig.savefig(out, dpi=180, facecolor=fig.get_facecolor())
    plt.close(fig)
    log.info("  -> %s", out)
    z = t.z.dropna()
    return {"parcels": len(t), "with_z": int(len(z)),
            "median_z": round(float(z.median()), 2) if len(z) else None,
            "share_below_minus1": round(float((z < -1).mean()), 3) if len(z) else None,
            "by_regime": {r: round(float(s.z.median()), 2)
                          for r, s in t.dropna(subset=["z"]).groupby("regime")},
            "file": str(out)}
