"""Cumulative precipitation and climatic water balance, day by day.

Key design decisions (vs. a naive monthly-cumsum plot):
  * daily resolution, so the current year is truncated at exactly the same
    day-of-year as the baseline -> no partial-month artefact;
  * the baseline is derived from the same station and same instrument network,
    not hard-coded and not taken from a different source;
  * the 30 individual baseline years are kept, so a percentile envelope shows
    whether the current year is outside the historical range, not merely
    "below average";
  * cumulative climatic water balance P - ET0 (Hargreaves), because
    precipitation alone understates a hot drought.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from drought.core.config import CFG

__all__ = ["BASELINE", "CURRENT_YEAR", "cumulative", "doy_matrix", "hargreaves_pet",
           "report_gaps", "stats"]

log = logging.getLogger(__name__)
BASELINE = tuple(CFG.precip.baseline)
CURRENT_YEAR = CFG.project.current_year


def doy_matrix(series: pd.Series, years: range) -> pd.DataFrame:
    """Wide frame: rows = calendar-aligned DOY 1..365, cols = year.

    29 February is dropped and later DOYs shifted, so a given row means the same
    calendar date in leap and non-leap years alike.
    """
    df = series.to_frame("v")
    df["year"] = df.index.year
    df["doy"] = df.index.dayofyear
    leap = df.index.is_leap_year
    after_feb29 = df["doy"] > 60
    df = df[~(leap & (df["doy"] == 60))]                       # drop 29 Feb
    df.loc[leap & after_feb29, "doy"] -= 1                     # realign the rest
    wide = df.pivot_table(index="doy", columns="year", values="v", aggfunc="first")
    return wide.reindex(index=range(1, 366), columns=list(years))


def cumulative(wide: pd.DataFrame) -> pd.DataFrame:
    """Cumulative sum per year; NaN days are treated as 0 but reported upstream."""
    return wide.fillna(0.0).cumsum()


def report_gaps(wide: pd.DataFrame, label: str, upto: int | None = None) -> None:
    sub = wide.iloc[:upto] if upto else wide
    gaps = sub.isna().sum()
    bad = gaps[gaps > 0]
    if not bad.empty:
        log.info("  [!] %s: missing days per year -> %s", label, bad.to_dict())
    else:
        log.info("  %s: no missing days", label)


# --------------------------------------------------------------------------- #
# Hargreaves reference evapotranspiration
# --------------------------------------------------------------------------- #
def hargreaves_pet(tmean, tmax, tmin, doy, lat_deg):
    """Reference ET0 [mm/day], Hargreaves-Samani (FAO-56 Ch. 4 fallback method)."""
    phi = np.radians(lat_deg)
    j = np.asarray(doy, dtype=float)
    dr = 1 + 0.033 * np.cos(2 * np.pi * j / 365.0)
    delta = 0.409 * np.sin(2 * np.pi * j / 365.0 - 1.39)
    ws = np.arccos(np.clip(-np.tan(phi) * np.tan(delta), -1.0, 1.0))
    ra_mj = (24 * 60 / np.pi) * 0.0820 * dr * (
        ws * np.sin(phi) * np.sin(delta) + np.cos(phi) * np.cos(delta) * np.sin(ws)
    )
    ra_mm = 0.408 * ra_mj
    dt = np.clip(np.asarray(tmax, float) - np.asarray(tmin, float), 0.0, None)
    return np.clip(0.0023 * ra_mm * (np.asarray(tmean, float) + 17.8) * np.sqrt(dt), 0.0, None)




def stats(cum_base: pd.DataFrame, cum_now: np.ndarray, n_days: int,
          signed: bool = False) -> tuple[float, float, float, int]:
    """Current value, baseline median, % of median and rank (1 = lowest) at DOY n_days.

    `signed=True` for quantities that can cross zero (P - ET0). A "% of median"
    is only meaningful for a strictly positive quantity: -84 against a median of
    +324 is NOT "-26% of normal", it is 408 mm below normal. Signed quantities
    therefore get no percentage (NaN).
    """
    now = float(cum_now[n_days - 1])
    med = float(cum_base.quantile(0.5, axis=1).iloc[n_days - 1])
    row = cum_base.iloc[n_days - 1].dropna()
    rank = int((row < now).sum()) + 1
    pct = float("nan") if (signed or med <= 0 or now < 0) else 100 * now / med
    return now, med, pct, rank
