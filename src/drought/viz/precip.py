"""The precipitation figure: cumulative P and P - ET0 against the climatology."""

from __future__ import annotations

import logging
import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from drought.analysis.precip import (
    BASELINE,
    CURRENT_YEAR,
    cumulative,
    doy_matrix,
    hargreaves_pet,
    report_gaps,
    stats,
)
from drought.core.config import CFG, figures_dir
from drought.sources.meteoswiss import P_PRECIP, P_TMAX, P_TMEAN, P_TMIN

__all__ = ["figure"]

log = logging.getLogger(__name__)
INCLUDE_WATER_BALANCE = CFG.precip.include_water_balance


MONTH_STARTS = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
MONTH_LABELS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]


def panel(ax, cum_base, cum_now, n_days, title, ylabel, unit="mm"):
    x = np.arange(1, 366)
    p10, p25, p50, p75, p90 = [cum_base.quantile(q, axis=1) for q in (.1, .25, .5, .75, .9)]

    ax.fill_between(x, p10, p90, color="#8b949e", alpha=0.20, linewidth=0,
                    label=f"{BASELINE[0]}-{BASELINE[1]} 10th-90th pct.")
    ax.fill_between(x, p25, p75, color="#8b949e", alpha=0.32, linewidth=0,
                    label=f"{BASELINE[0]}-{BASELINE[1]} 25th-75th pct.")
    ax.plot(x, p50, color="white", ls="--", lw=2.0,
            label=f"{BASELINE[0]}-{BASELINE[1]} median")
    ax.plot(x[:n_days], cum_now[:n_days], color="#ff3b30", lw=2.6,
            label=f"{CURRENT_YEAR} (to DOY {n_days})")
    ax.plot([n_days], [cum_now[n_days - 1]], "o", color="#ff3b30", ms=5)

    ax.set_title("\n".join(textwrap.wrap(title, width=40)),
                 fontsize=12, fontweight="bold", color="white", pad=10)
    ax.set_ylabel(f"{ylabel} [{unit}]", fontsize=10, color="white")
    ax.set_xlim(1, 365)
    ax.set_xticks(MONTH_STARTS)
    ax.set_xticklabels(MONTH_LABELS, color="white")
    ax.tick_params(colors="white", labelsize=9, labelbottom=True)
    ax.grid(True, ls=":", alpha=0.25, color="gray")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#484f58")
    ax.legend(frameon=True, facecolor="#161b22", edgecolor="none",
              fontsize=8.5, labelcolor="white", loc="upper left")
    ax.set_box_aspect(1)   # square plotting area, independent of label sizes


def annotate(ax, cum_base, cum_now, n_days, unit="mm", signed=False):
    """Annotate with the current value, the baseline median, and the rank.

    `signed=True` for quantities that can cross zero (P - ET0). A "% of median"
    is only meaningful for a strictly positive quantity: -84 against a median of
    +324 is NOT "-26% of normal", it is 408 mm below normal. Reporting a ratio
    there invites a serious misreading, so signed quantities get an absolute
    difference instead.
    """
    now, med, pct, rank = stats(cum_base, cum_now, n_days, signed)
    comp = f"{now - med:+.0f} {unit} vs. median" if np.isnan(pct) else f"{pct:.0f}% of median"
    n = int(cum_base.iloc[n_days - 1].notna().sum()) + 1
    txt = (f"{now:.0f} {unit}   |   median {med:.0f} {unit}   "
           f"({comp})   |   rank {rank}/{n}")
    ax.text(0.985, 0.04, txt, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=9, color="#c9d1d9",
            bbox=dict(facecolor="#161b22", edgecolor="#30363d", boxstyle="round,pad=0.4"))
    return now, med, pct, rank


def figure(daily: pd.DataFrame, abbr: str, lat: float,
           station_label: str | None = None) -> dict:
    """Cumulative precipitation (and P - ET0) figure for the current year.

    Args:
        daily: The station's daily series (sources.meteoswiss.load_local).
        abbr: Station abbreviation, for the file name.
        lat: Station latitude, for ET0.
        station_label: Name on the figure (default precip.station_name).

    Returns:
        dict with end_doy, precip and water-balance totals, medians and ranks.
    """
    label = station_label or CFG.precip.station_name
    log.info("  %s daily records, %s -> %s", f"{len(daily):,}", daily.index.min().date(),
             daily.index.max().date())

    base_years = range(BASELINE[0], BASELINE[1] + 1)

    # --- precipitation -----------------------------------------------------
    p_base = doy_matrix(daily[P_PRECIP], base_years)
    p_now = doy_matrix(daily[P_PRECIP], range(CURRENT_YEAR, CURRENT_YEAR + 1))
    n_days = int(p_now[CURRENT_YEAR].last_valid_index() or 0)
    if n_days == 0:
        raise ValueError(f"No {CURRENT_YEAR} data at station {abbr!r}.")
    log.info("  current year complete through aligned DOY %d", n_days)
    report_gaps(p_base, "baseline precipitation")
    report_gaps(p_now, f"{CURRENT_YEAR} precipitation", upto=n_days)

    cum_p_base = cumulative(p_base)
    cum_p_now = cumulative(p_now)[CURRENT_YEAR].to_numpy()

    # --- optional water balance -------------------------------------------
    have_wb = INCLUDE_WATER_BALANCE and all(
        c in daily.columns for c in (P_TMEAN, P_TMAX, P_TMIN)
    )
    if INCLUDE_WATER_BALANCE and not have_wb:
        log.warning("  [!] temperature parameters unavailable - skipping water-balance panel")

    if have_wb:
        pet = pd.Series(
            hargreaves_pet(daily[P_TMEAN], daily[P_TMAX], daily[P_TMIN],
                           daily.index.dayofyear, lat),
            index=daily.index,
        )
        wb = daily[P_PRECIP] - pet
        wb_base = doy_matrix(wb, base_years)
        wb_now = doy_matrix(wb, range(CURRENT_YEAR, CURRENT_YEAR + 1))
        cum_wb_base = cumulative(wb_base)
        cum_wb_now = cumulative(wb_now)[CURRENT_YEAR].to_numpy()

    # --- plot ---------------------------------------------------------------
    plt.style.use("dark_background")
    nrows = 2 if have_wb else 1
    fig, axes = plt.subplots(nrows, 1, figsize=(6.8, 7.6 * nrows), dpi=200,
                             sharex=True, constrained_layout=True)
    axes = np.atleast_1d(axes)

    panel(axes[0], cum_p_base, cum_p_now, n_days,
          f"{label} - cumulative precipitation", "Cumulative precipitation")
    st = annotate(axes[0], cum_p_base, cum_p_now, n_days)
    log.info("\nPrecipitation to DOY %d: %.0f mm (%.0f%% of the %d-%d median), rank %d of %d "
             "IN THE COMPARISON WINDOW (not the full station record)", n_days, st[0], st[2],
             BASELINE[0], BASELINE[1], st[3], BASELINE[1] - BASELINE[0] + 2)

    if have_wb:
        panel(axes[1], cum_wb_base, cum_wb_now, n_days,
              f"{label} - cumulative water balance P-ET0 (Hargreaves)",
              "Cumulative P - ET0")
        s2 = annotate(axes[1], cum_wb_base, cum_wb_now, n_days, signed=True)
        log.info("Water balance to DOY %d: %.0f mm (median %.0f mm), rank %d",
                 n_days, s2[0], s2[1], s2[3])

    # The credit is the x label (see viz.condition.figure).
    axes[-1].set_xlabel("Data: MeteoSchweiz Open Data (SwissMetNet), CC-BY",
                        loc="left", fontsize=7.5, color="#6e7681")
    out = figures_dir("station") / f"precip_{abbr}_{CURRENT_YEAR}.png"
    fig.savefig(out, dpi=200, facecolor=fig.get_facecolor())
    log.info("  -> %s", out)
    plt.close(fig)
    res = {"end_doy": n_days, "precip_mm": st[0], "precip_median_mm": st[1],
           "precip_pct_of_median": st[2], "precip_rank": st[3],
           "n_years": BASELINE[1] - BASELINE[0] + 2, "file": str(out)}
    if have_wb:
        res.update(wb_mm=s2[0], wb_median_mm=s2[1], wb_rank=s2[3])
    return res



