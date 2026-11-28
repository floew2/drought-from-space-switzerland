"""The VHI cross-check figure: swisstopo's VHI on our pure grassland cells."""

from __future__ import annotations

import logging

import matplotlib.pyplot as plt
import pandas as pd

from drought.core.config import CFG, figures_dir, series_dir
from drought.viz import credit

__all__ = ["figure", "load"]

log = logging.getLogger(__name__)


def load(site: str, what: str = "series") -> pd.DataFrame:
    """A site's VHI table ("series" or "parcels"), verified."""
    from drought.core.provenance import verify

    f = series_dir(site) / f"vhi_{what}.csv"
    if not f.exists():
        raise FileNotFoundError(f"no {f}; run `drought vhi --aoi {site}`")
    verify(f, "vhi")
    return pd.read_csv(f, dtype={"parcel_id": str})


def figure(sites: list[str]) -> dict:
    """One panel per site; returns the season mean per site and system."""
    cur = int(CFG.project.current_year)
    fig, axes = plt.subplots(1, len(sites), figsize=(4.2 * len(sites), 3.4),
                             sharey=True, squeeze=False)
    out: dict = {}
    for ax, site in zip(axes[0], sites):
        ax.axhspan(40, 50, color="#cbffca", alpha=0.6, lw=0, label="normal (40-50)")
        ax.axhspan(0, 30, color="#f5cd85", alpha=0.35, lw=0, label="stressed (< 30)")
        try:
            d = load(site)
        except FileNotFoundError as exc:
            log.warning("  [skipped] %s", exc)
            continue
        d["date"] = pd.to_datetime(d.date)
        for r, g in d.groupby("regime"):
            ax.plot(g.date, g.vhi_mean, lw=1.4, label=r)
            out[f"{site}/{r}"] = {"vhi_mean": round(float(g.vhi_mean.mean()), 1),
                                  "days": int(len(g))}
        ax.set_title(CFG.sites.get(site).label, fontsize=10)
        ax.set_ylim(0, 100)
        ax.tick_params(axis="x", labelrotation=30, labelsize=8)
    axes[0][0].set_ylabel("VHI (swissEO)")
    axes[0][-1].legend(fontsize=7, loc="lower right")
    fig.suptitle(f"Vegetation Health Index on pure grassland cells, {cur}", fontsize=11)
    fig.supxlabel("VHI 0-100, low is stressed.", fontsize=7, color="#6e7681")
    fig.tight_layout()
    f = figures_dir() / f"vhi_{cur}.png"
    credit.save(fig, f, "vhi", "parcels", dpi=150)
    plt.close(fig)
    log.info("  -> %s", f)
    return out
