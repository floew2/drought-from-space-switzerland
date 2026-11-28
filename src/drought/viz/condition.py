"""Condition figures: the season trajectory and the season integral per site
and management system, in the same visual grammar as the precipitation figure:
percentile envelope, white dashed median, red current year, square panels.

Panel A: seasonal index trajectory. Shows WHEN 2026 diverged, which a single
         anomaly number cannot.
Panel B: cumulative integral of the index. The analogue of cumulative
         precipitation; integrated index is an established proxy for season
         biomass production.
"""

from __future__ import annotations

import logging
import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from drought.analysis import condition as C
from drought.analysis.condition import (
    BASELINE,
    CURRENT_YEAR,
    DOY_MAX,
    DOY_MIN,
    INTEGRATION_WINDOW,
    MONTH_LABELS,
    MONTH_STARTS,
    REGIME_LABELS,
    SMOOTH_DAYS,
)
from drought.core.config import CFG, figures_dir
from drought.viz import credit

__all__ = ["figure", "plot"]

log = logging.getLogger(__name__)


def smooth(v: np.ndarray, days: int) -> np.ndarray:
    if not days:
        return v
    return (pd.Series(v).rolling(days, center=True, min_periods=max(2, days // 3))
            .median().to_numpy())


def panel(ax, base: pd.DataFrame, cur: np.ndarray, ylabel: str, first: bool,
          raw: tuple | None = None):
    x = base.index.to_numpy()
    p10, p25, p50, p75, p90 = [base.quantile(q, axis=1) for q in (.1, .25, .5, .75, .9)]

    lo, hi = BASELINE
    ax.fill_between(x, p10, p90, color="#8b949e", alpha=0.20, lw=0,
                    label=f"{lo}-{hi} 10th-90th pct.")
    ax.fill_between(x, p25, p75, color="#8b949e", alpha=0.32, lw=0,
                    label=f"{lo}-{hi} 25th-75th pct.")
    # The smoothing note belongs only where smoothing was applied: the
    # trajectory panel passes `raw`, the cumulative panel does not.
    sm = f" ({SMOOTH_DAYS}-day rolling median)" if (SMOOTH_DAYS and raw is not None) else ""
    ax.plot(x, p50, color="white", ls="--", lw=2.0,
            label=f"{lo}-{hi} median{sm}")
    ax.plot(x, cur, color="#ff3b30", lw=2.4, label=f"{CURRENT_YEAR}{sm}")
    if raw is not None:
        # Actual observations, so smoothing never hides what was measured.
        ax.plot(raw[0], raw[1], "o", color="#ff3b30", ms=2.6, alpha=0.55,
                lw=0, label=f"{CURRENT_YEAR} observations (raw)")

    ax.set_box_aspect(1)
    ax.set_ylabel(ylabel, fontsize=10, color="white")
    ax.set_xlim(DOY_MIN, DOY_MAX)
    ax.set_xticks(MONTH_STARTS)
    ax.set_xticklabels(MONTH_LABELS, color="white")
    ax.tick_params(colors="white", labelsize=9, labelbottom=True)
    ax.grid(True, ls=":", alpha=0.25, color="gray")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#484f58")
    if first:
        ax.legend(frameon=True, facecolor="#161b22", edgecolor="none",
                  fontsize=8, labelcolor="white", loc="lower left")


def figure(site: str, regime: str, sub: pd.DataFrame):
    """The two-panel figure for one site and system; returns its numbers."""
    p = C.integral_parts(sub)
    if p is None:
        log.info("    [!] insufficient data for %s/%s — skipped", site, regime)
        return None
    lo, hi = INTEGRATION_WINDOW
    base, base_years, days = p["base"], p["base_years"], p["days"]
    cum_base, cum_cur, cur = p["cum_base"], p["cum_cur"], p["cur"]
    n_days_used = int(days.sum())
    log.info("    window DOY %d-%d, %d baseline years, %d runs to DOY %d", lo, hi,
             len(base_years), CURRENT_YEAR,
             int(p["wide"].index[np.where(~np.isnan(cur))[0][-1]]))
    log.info("    integral over %d days shared by all %d years", n_days_used,
             len(base_years) + 1)
    if n_days_used < 40:
        log.info("    [!] few shared days — the integral is thin, prefer the "
                 "trajectory panel and the mean anomaly")

    plt.style.use("dark_background")
    fig, axes = plt.subplots(2, 1, figsize=(6.4, 12.4), dpi=180,
                             constrained_layout=True)

    title = f"{CFG.sites.get(site).label} - {REGIME_LABELS.get(regime, regime)}"
    obs_doy = sub.loc[sub["year"] == CURRENT_YEAR, "doy"].to_numpy()
    obs_val = sub.loc[sub["year"] == CURRENT_YEAR, C.INDEX].to_numpy()
    base_s = base.apply(lambda c: smooth(c.to_numpy(), SMOOTH_DAYS))
    base_s.index = base.index
    panel(axes[0], base_s, smooth(cur, SMOOTH_DAYS), C.INDEX, first=True,
          raw=(obs_doy, obs_val))
    axes[0].set_title("\n".join(textwrap.wrap(title, 34)), fontsize=12,
                      fontweight="bold", color="white", pad=10)
    panel(axes[1], cum_base, cum_cur,
          f"Cumulative {C.INDEX} (season integral)", first=False)

    pct, rank, upto = p["pct"], p["rank"], p["upto"]
    # Shade everything not contributing to the integral.
    excl = ~days
    axes[1].fill_between(base.index.to_numpy(), 0, 1, where=excl,
                         transform=axes[1].get_xaxis_transform(),
                         color="#000000", alpha=0.45, lw=0)
    axes[1].text(0.97, 0.06,
                 f"integral to DOY {upto}: {pct:+.0f}% vs. median   |   "
                 f"rank {rank}/{p['n_years']}\n{n_days_used} shared days "
                 f"in DOY {lo}-{hi}",
                 transform=axes[1].transAxes, ha="right", fontsize=9,
                 color="#c9d1d9",
                 bbox=dict(facecolor="#161b22", edgecolor="#30363d",
                           boxstyle="round,pad=0.4"))

    out = credit.save(fig, figures_dir(site) / f"vi_trajectory_{C.INDEX}_{regime}.png",
                      "s2", "parcels")
    log.info("  -> %s", out)
    plt.close(fig)

    n_obs = {int(y): len(o) for y, o in p["observed"].items()}
    return dict(site=site, regime=regime, integral_pct=pct, rank=rank,
                upto=upto, window=(lo, hi), n_days=n_days_used,
                n_baseline_years=len(base_years), obs=n_obs, file=out)




def plot(sites=None, index=None, regimes=None, verbose=True) -> pd.DataFrame | None:
    """Trajectory and season-integral figures, one per site and system.

    The trajectory panel is the trustworthy one: it shows WHEN the divergence
    happened and depends on no interpolation choices. The integral is sensitive
    to how gaps are handled, which is why it is the second panel and not the
    headline.

    Returns:
        One row per site and system with the integral's numbers, or None.
    """
    if index:
        C.INDEX = index
    out = []
    for s_ in sites or list(CFG.condition.sites):
        df = C.load(s_, verbose=verbose)
        for regime, g in df.groupby("regime"):
            if regimes and regime not in regimes:
                continue
            log.info("\n%s / %s: %d obs", s_, regime, len(g))
            r = figure(s_, regime, g)
            if r:
                out.append(r)
    return pd.DataFrame(out) if out else None
