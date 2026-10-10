"""Synthetic MeteoSwiss station series: daily precipitation and temperature.

Written in the MeteoSwiss Open Data CSV layout (semicolon-separated,
`reference_timestamp` as dd.mm.yyyy HH:MM, parameters rre150d0, tre200d0,
tre200dx, tre200dn) so `precip` reads it through the same code path as the
real files.

Calibrated to the real result (targets.PRECIP): the 1991-2020 median
cumulative precipitation to 10 August is 826 mm and 2026 reaches 501 mm, the
driest of 31; the Hargreaves water balance P - ET0 has a median of +324 mm and
2026 ends at -84 mm. Precipitation is scaled to hit the totals exactly; ET0 is
hit by solving for a temperature-range factor (baseline) and a summer heat
anomaly (2026, on top of a 30 % wider daily range under clear skies).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import optimize

from . import targets as T

ABBR = "syn"
LAT = 47.425                    # St. Gallen
ELEV = 776
BASELINE = (1991, 2020)
CURRENT = 2026
# Mean monthly precipitation, mm — the climatological shape at a north-Alpine
# foothill station: wet summers, drier winters. Rescaled during calibration.
MONTHLY_MM = [70, 65, 80, 90, 125, 150, 170, 165, 115, 90, 85, 80]


def _monthly_mean_daily(days: pd.DatetimeIndex) -> np.ndarray:
    m = days.month.to_numpy() - 1
    return np.array(MONTHLY_MM, float)[m] / days.days_in_month.to_numpy()


def _precip(rng, days: pd.DatetimeIndex) -> np.ndarray:
    """Markov-chain occurrence, gamma amounts."""
    mean = _monthly_mean_daily(days)
    p_wet_wet, p_wet_dry = 0.62, 0.32
    pw = p_wet_dry / (1 - p_wet_wet + p_wet_dry)       # stationary wet share
    wet = np.zeros(len(days), bool)
    for i in range(1, len(days)):
        wet[i] = rng.uniform() < (p_wet_wet if wet[i - 1] else p_wet_dry)
    shape = 0.75
    amt = rng.gamma(shape, mean / pw / shape)
    return np.where(wet, amt, 0.0)


def _temperature(rng, days: pd.DatetimeIndex, dtr_scale: float, heat: np.ndarray,
                 sunny: np.ndarray):
    doy = days.dayofyear.to_numpy()
    seasonal = 8.0 + 9.5 * np.sin(2 * np.pi * (doy - 110) / 365.25)
    anom = np.zeros(len(days))
    for i in range(1, len(days)):
        anom[i] = 0.8 * anom[i - 1] + rng.normal(0, 1.6)
    tmean = seasonal + anom + heat
    dtr = dtr_scale * (7.5 + 2.0 * np.sin(2 * np.pi * (doy - 100) / 365.25)) \
        * np.exp(rng.normal(0, 0.15, len(days))) * (1.0 + 0.3 * sunny)
    return tmean, tmean + dtr / 2, tmean - dtr / 2


def _cum_to(series: pd.Series, year: int, end_doy: int) -> float:
    s = series[(series.index.year == year) & (series.index.dayofyear <= end_doy)]
    return float(s.sum())


def generate(rng: np.random.Generator) -> pd.DataFrame:
    """Daily station series 1991-01-01 .. 2026-08-10, calibrated.

    Returns:
        DataFrame indexed by date with rre150d0, tre200d0, tre200dx, tre200dn.
    """
    from ..precip import hargreaves_pet

    tgt = T.PRECIP
    end = tgt["end_doy"]
    days = pd.date_range(f"{BASELINE[0]}-01-01", f"{CURRENT}-12-31", freq="D")
    p = pd.Series(_precip(rng, days), index=days)
    base_years = range(BASELINE[0], BASELINE[1] + 1)

    # Precipitation: scale the baseline to the median, then 2026 to its total.
    med = np.median([_cum_to(p, y, end) for y in base_years])
    p *= tgt["p_median_mm"] / med
    in26 = (p.index.year == CURRENT) & (p.index.dayofyear <= end)
    p[in26] *= tgt["p_2026_mm"] / _cum_to(p, CURRENT, end)
    # Spring near normal, summer dry: move rain out of June-August 2026 into
    # its first months while keeping the 2026 total — then re-check the total.
    summer = in26 & (p.index.month >= 6)
    early = in26 & (p.index.month <= 4)
    moved = 0.35 * p[summer].sum()
    p[summer] *= 0.65
    p[early] += moved * p[early] / p[early].sum()
    p[(p.index.year == CURRENT) & (p.index.dayofyear > end)] = np.nan

    noise_rng = np.random.default_rng(int(rng.integers(1 << 31)))
    state = noise_rng.bit_generator.state
    summer26 = ((days.year == CURRENT) & (days.month >= 5) & (days.month <= 8)).astype(float)

    def series(dtr_scale, heat_c):
        noise_rng.bit_generator.state = state           # same weather, rescaled
        # 2026 summer: hotter, and sunnier — clear skies widen the daily range.
        tm, tx, tn = _temperature(noise_rng, days, dtr_scale, heat_c * summer26, summer26)
        et0 = pd.Series(hargreaves_pet(tm, tx, tn, days.dayofyear, LAT), index=days)
        return tm, tx, tn, et0

    def wb_median(dtr_scale):
        # the 2026 terms do not touch baseline years
        et0 = series(dtr_scale, 0.0)[3]
        return np.median([_cum_to(p, y, end) - _cum_to(et0, y, end) for y in base_years])

    dtr_scale = optimize.brentq(lambda s: wb_median(s) - tgt["wb_median_mm"], 0.3, 3.0)

    def wb_2026(heat):
        et0 = series(dtr_scale, heat)[3]
        return _cum_to(p, CURRENT, end) - _cum_to(et0, CURRENT, end)

    heat = optimize.brentq(lambda h: wb_2026(h) - tgt["wb_2026_mm"], -5.0, 12.0)
    tm, tx, tn, _ = series(dtr_scale, heat)
    df = pd.DataFrame({"rre150d0": p.round(1), "tre200d0": np.round(tm, 1),
                       "tre200dx": np.round(tx, 1), "tre200dn": np.round(tn, 1)},
                      index=days)
    df = df[~((df.index.year == CURRENT) & (df.index.dayofyear > end))]
    df.attrs.update(dtr_scale=float(dtr_scale), heat_2026_c=float(heat))
    return df


def write(df: pd.DataFrame, folder) -> list:
    """Write the series as MeteoSwiss-style 'historical' and 'recent' files.

    Rounding to 0.1 mm / 0.1 degC happens here, as in the real files, so the
    totals read back can differ from the calibration targets by a fraction of
    a millimetre.
    """
    from pathlib import Path

    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    out = df.copy()
    out.insert(0, "reference_timestamp", out.index.strftime("%d.%m.%Y 00:00"))
    out.insert(0, "station_abbr", ABBR.upper())
    files = []
    for name, sel in (("historical", out.index.year < 2020),
                      ("recent", out.index.year >= 2020)):
        f = folder / f"ogd-smn_{ABBR}_d_{name}.csv"
        out[sel].to_csv(f, sep=";", index=False)
        files.append(f)
    meta = pd.DataFrame([{"station_abbr": ABBR.upper(), "station_name": "Synthetic station",
                          "station_coordinates_wgs84_lat": LAT,
                          "station_height_masl": ELEV,
                          "station_data_since": f"01.01.{BASELINE[0]}"}])
    meta.to_csv(folder / "ogd-smn_meta_stations.csv", sep=";", index=False)
    return files
