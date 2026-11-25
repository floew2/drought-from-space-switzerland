"""Grassland condition: cleaning the site series and the season integral.

Input: <paths.series>/<site>/condition.csv (products/condition.py), columns
site, regime, year, date, doy, NDVI, kNDVI, NDMI, NDRE, n_px.

`prepare` removes residual haze and partial-coverage dates; `season_integral`
compares the current year's cumulative index with the baseline years over the
days every year covers. The figures are in viz/condition.py.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from drought.core.config import CFG, series_dir

log = logging.getLogger(__name__)

# Which index to plot. The CSV carries all four; switching is free.
#   kNDVI  recommended default — same information as NDVI but without the
#          saturation that compresses dense Swiss swards.
#   NDMI   canopy water; earlier signal, but snow inflates it (alpine caution).
#   NDRE   red-edge chlorophyll; Sentinel-2's real advantage over MODIS.
#   NDVI   familiar, weakest discriminator here.
INDEX = CFG.detect.index
ALL_INDICES = ["NDVI", "kNDVI", "NDMI", "NDRE"]
SMOOTH_DAYS = CFG.analysis.smooth_days              # rolling-median window on the daily grid, 0 = off.
                              # A single hazy or partly-covered date should not
                              # be readable as a real drop; the raw observations
                              # stay visible as dots so nothing is hidden.

# --- residual-haze filter --------------------------------------------------- #
# A cloud mask lets some thin haze and cirrus through. Its
# fingerprint is unmistakable and wavelength-ordered: scattering inflates the
# short-wavelength bands, so RED is hit hardest, RED-EDGE less, and SWIR almost
# not at all. Rheintal 2026 lost 0.197 NDVI, 0.086 NDRE and only 0.028 NDMI
# between DOY 120 and 122 — two days, full site coverage. No vegetation change
# does that; mowing would not spare NDMI either.
#
# So: flag dates where NDVI sits far below its own local level while NDMI does
# not follow. Using the DISAGREEMENT between indices means no absolute
# threshold has to be guessed.
HAZE_FILTER = True
HAZE_WINDOW = int(CFG.clean.haze_window)  # neighbouring observations for the local median
HAZE_NDVI_Z = float(CFG.clean.haze_ndvi_z) # how far NDVI must fall below its local level (robust z)
HAZE_GAP = float(CFG.clean.haze_gap)       # how much less NDMI fell, in the same units
CURRENT_YEAR = CFG.project.current_year
BASELINE = (CFG.project.baseline_years[0], CFG.project.baseline_years[-1])

DOY_MIN, DOY_MAX = (int(v) for v in CFG.condition.doy)

# FIXED integration window for the season integral.
# The earlier data-driven "common to all years" window was a mistake: it took
# the intersection across years, so a single cloudy September in one baseline
# year truncated the window for everyone — closing the integral at mid-August,
# right where the 2026 deficit begins, while still accumulating two months of an
# unusually strong early spring beforehand. That produced positive integrals
# alongside negative growing-season means.
# 1 May - 30 Sep covers the growing season without the spring shoulder.
INTEGRATION_WINDOW = tuple(CFG.analysis.integration_window)
MIN_WINDOW_COVERAGE = CFG.analysis.min_window_coverage     # a year needs this fraction of the window covered
                              # by interpolated (non-gap) values to be counted
MIN_COVERAGE = CFG.analysis.min_coverage            # fraction of the site's mask that must be clear.
                              # Raised from 0.5: a date where only half the site
                              # is cloud-free samples a DIFFERENT set of parcels
                              # than a full date, and if cut timing or elevation
                              # vary spatially the mean shifts sharply. That is
                              # the leading suspect for large day-to-day jumps.
                              # An ABSOLUTE pixel threshold does not work: the
                              # mown and grazed masks differ in size by ~7x, so
                              # one number cannot serve both.
MIN_OBS_PER_YEAR = CFG.analysis.min_obs_per_year          # below this a year cannot support a curve
MAX_GAP_DAYS = CFG.analysis.max_gap_days             # do not interpolate across gaps longer than this

REGIME_LABELS = {"mown": "mown grassland", "grazed": "grazed pasture",
                 "intensive": "intensive meadows", "extensive": "extensive meadows",
                 "permanent": "permanent grassland"}

# Management systems are kept SEPARATE (config: analysis.merge_regimes). The
# trajectory export carries intensive and extensive meadows; pooling them
# averages two cutting regimes into a curve that belongs to neither. The pooled
# "permanent grassland" series is still available, area-weighted by n_px.
MERGE_REGIMES = bool(getattr(CFG.analysis, "merge_regimes", False))

MONTH_STARTS = [60, 91, 121, 152, 182, 213, 244, 274]
MONTH_LABELS = ["Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct"]


# --------------------------------------------------------------------------- #
# Loading and gridding
# --------------------------------------------------------------------------- #
def path(site: str):
    return series_dir(site) / "condition.csv"


def load(site: str, verbose: bool = True) -> pd.DataFrame:
    """A site's condition series, verified against the current settings, cleaned."""
    from drought.core.provenance import verify

    f = path(site)
    if not f.exists():
        raise FileNotFoundError(f"no {f}; run `drought condition --aoi {site}`")
    verify(f, "condition")
    if verbose:
        log.info("reading %s", f)
    return prepare(pd.read_csv(f), verbose=verbose)


def prepare(df: pd.DataFrame, verbose: bool = True) -> pd.DataFrame:
    """Clean a trajectory export: tile duplicates, haze, partial coverage.

    Args:
        df: Rows in the export schema (site, regime, year, date, doy, the
            four indices, n_px), as read from CSV.
        verbose: Print what each filter removed.

    Returns:
        Cleaned frame with `ix_mean` (the selected index) and `cov`.
    """
    df = df.copy()
    df.columns = [c.strip().lower() for c in df.columns]
    for c in ("year", "doy", "n_px"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    present = [c for c in ALL_INDICES if c.lower() in df.columns]
    df = df.rename(columns={c.lower(): c for c in ALL_INDICES})
    if INDEX not in present:
        raise KeyError(f"{INDEX!r} not in CSV; found {present}")
    for c in present:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=present + ["doy", "year"])

    def wmean(g, how):
        """n_px-weighted mean of every index; n_px combined per `how`."""
        w = g["n_px"].to_numpy(float)
        out = {c: float(np.average(g[c], weights=w))
               for c in ALL_INDICES if c in g.columns}
        out["n_px"] = float(w.max() if how == "max" else w.sum())
        return pd.Series(out)

    # ---- stage A: collapse MGRS tile duplicates, within a regime -----------
    # Sentinel-2 tiles overlap, so a site on a tile boundary produces two rows
    # for one date covering LARGELY THE SAME GROUND. The index is a weighted
    # mean (a 150-pixel sliver must not count as much as a 6200-pixel scene),
    # but n_px takes the MAX, not the sum: summing double-counts the overlap,
    # and a single date with both tiles clear would then set an inflated maximum
    # that makes every ordinary date look half-covered.
    df = (df.groupby(["site", "regime", "year", "doy"], as_index=False)
            .apply(lambda g: wmean(g, "max"), include_groups=False)
            .reset_index(drop=True))

    # ---- stage B: pool regimes into one permanent-grassland series ---------
    # Summing IS correct here: mown and grazed are disjoint areas on the ground.
    if MERGE_REGIMES:
        df = (df.groupby(["site", "year", "doy"], as_index=False)
                .apply(lambda g: wmean(g, "sum"), include_groups=False)
                .reset_index(drop=True))
        df["regime"] = "permanent"

    if HAZE_FILTER and {"NDVI", "NDMI"} <= set(df.columns):
        def mark(g):
            g = g.sort_values("doy").copy()
            z = {}
            for ix in ("NDVI", "NDMI"):
                res = g[ix] - g[ix].rolling(HAZE_WINDOW, center=True,
                                            min_periods=2).median()
                mad = (res - res.median()).abs().median() * 1.4826
                z[ix] = res / max(mad, 1e-3)
            g["haze"] = (z["NDVI"] < HAZE_NDVI_Z) & \
                        ((z["NDMI"] - z["NDVI"]) > HAZE_GAP)
            return g

        df = (df.groupby(["site", "regime", "year"], group_keys=False)
                .apply(mark, include_groups=False)
                .join(df[["site", "regime", "year"]]))
        n_haze = int(df["haze"].sum())
        by_year = df[df["haze"]].groupby("year").size()
        if verbose:
            log.info("  haze filter removed %d of %d dates (%.1f%%)", n_haze, len(df),
                     100 * n_haze / max(len(df), 1))
        if verbose and len(by_year):
            log.info("    by year: %s  <- roughly even means the filter is not biasing a year",
                     ", ".join(f"{int(y)}:{n}" for y, n in by_year.items()))
        df = df[~df["haze"]].drop(columns="haze")

    df["ix_mean"] = df[INDEX]

    # Coverage relative to the best observation for that site and regime.
    full = df.groupby(["site", "regime"])["n_px"].transform("max")
    df["cov"] = df["n_px"] / full

    before = len(df)
    df = df[df["cov"] >= MIN_COVERAGE]
    if verbose:
        log.info("  dropped %d of %d dates below %.0f%% mask coverage", before - len(df),
                 before, 100 * MIN_COVERAGE)
    return df


def to_grid(sub: pd.DataFrame, verbose: bool = True) -> tuple[pd.DataFrame, dict]:
    """Interpolate each year onto a common DOY axis.

    Cloud means every year has observations on different dates, so the years
    cannot be compared until they share an axis. Gaps longer than MAX_GAP_DAYS
    are left as NaN rather than bridged: a 6-week hole spanning a cut would
    otherwise be silently invented, and invented data in the 2026 curve is
    exactly what a sceptical reader will go looking for.
    """
    doy = np.arange(DOY_MIN, DOY_MAX + 1)
    grid, observed = {}, {}
    for year, g in sub.groupby("year"):
        g = g.sort_values("doy")
        if len(g) < MIN_OBS_PER_YEAR:
            if verbose:
                log.info("    [!] %d: only %d clear dates — skipped", int(year), len(g))
            continue
        vals = np.interp(doy, g["doy"], g["ix_mean"], left=np.nan, right=np.nan)
        # blank out long gaps
        obs = g["doy"].to_numpy()
        for a, b in zip(obs[:-1], obs[1:]):
            if b - a > MAX_GAP_DAYS:
                vals[(doy > a) & (doy < b)] = np.nan
        vals[doy < obs[0]] = np.nan
        vals[doy > obs[-1]] = np.nan
        grid[int(year)] = vals
        observed[int(year)] = obs
    return pd.DataFrame(grid, index=doy), observed


def window_coverage(wide: pd.DataFrame, lo: int, hi: int) -> pd.Series:
    """Fraction of the fixed window each year actually covers (non-gap)."""
    sel = wide.loc[lo:hi]
    return sel.notna().sum() / len(sel)


def common_days(wide: pd.DataFrame, lo: int, hi: int) -> np.ndarray:
    """Boolean mask of days inside the window that EVERY year covers.

    Without this the integral is partly a measure of cloud luck: treating a gap
    as zero growth means a year with a long gap accumulates less, so a well
    observed year beats a poorly observed one regardless of vegetation. That
    produced the contradiction of a curve 15% below the median ranking only
    4th of 8.
    """
    inwin = (wide.index >= lo) & (wide.index <= hi)
    return inwin & wide.notna().all(axis=1).to_numpy()


def cumulative(wide: pd.DataFrame, days: np.ndarray) -> pd.DataFrame:
    """Season integral over exactly the days in `days`, for every year alike."""
    out = wide.copy()
    out.loc[~days] = 0.0
    return out.fillna(0.0).cumsum()


def integral_parts(sub: pd.DataFrame, verbose: bool = True) -> dict | None:
    """Everything the season integral needs, shared by the number and the figure."""
    wide, observed = to_grid(sub, verbose=verbose)
    base_years = [y for y in wide.columns if BASELINE[0] <= y <= BASELINE[1]]
    if CURRENT_YEAR not in wide.columns or len(base_years) < 3:
        if verbose:
            log.info("    [!] insufficient data — skipped")
        return None

    lo, hi = INTEGRATION_WINDOW
    cov = window_coverage(wide[base_years + [CURRENT_YEAR]], lo, hi)
    # The current year is partial by definition — the season is not over. That
    # is not a data-quality problem and must not be reported as one; the
    # comparison simply ends where 2026 ends.
    thin = cov[base_years][cov[base_years] < MIN_WINDOW_COVERAGE]
    if verbose and not thin.empty:
        log.info("    [!] gappy baseline years dropped: %s",
                 ", ".join(f"{int(y)}:{c:.0%}" for y, c in thin.items()))
    base_years = [y for y in base_years if cov[y] >= MIN_WINDOW_COVERAGE]
    if len(base_years) < 3:
        if verbose:
            log.info("    [!] fewer than 3 usable baseline years — skipped")
        return None
    if verbose and len(base_years) < 5:
        log.info("    [!] only %d baseline years survive — the percentile envelope is weak "
                 "here, treat the rank with caution", len(base_years))
    base = wide[base_years]
    cur = wide[CURRENT_YEAR].to_numpy()

    # Days covered by every year, current included, so all years integrate over
    # an identical set of dates.
    days = common_days(wide[base_years + [CURRENT_YEAR]], lo, hi)
    cum_base = cumulative(base, days)
    cum_cur = cumulative(wide[[CURRENT_YEAR]], days)[CURRENT_YEAR].to_numpy()

    # Stop the current-year line at its last real observation. Otherwise the
    # cumulative curve runs flat to the right edge, which reads as "growth
    # stopped" rather than "data ends here". Same convention as the
    # precipitation figure.
    valid_cur = np.where(~np.isnan(cur))[0]
    last_obs = valid_cur[-1] if len(valid_cur) else len(cur) - 1
    cum_cur = cum_cur.copy()
    cum_cur[last_obs + 1:] = np.nan

    # Compare at the last date the CURRENT year actually has data for, not at
    # the window end. Clipping cum_cur to the last observation left NaN beyond
    # it, and reading the headline number at DOY 273 therefore returned NaN
    # with a meaningless rank of 1. Every year is summed over the same days,
    # ending where 2026 stops — a season-to-date integral, not a full season.
    win_end = int(np.where(base.index.to_numpy() == hi)[0][0])
    last = min(win_end, last_obs)
    now = float(cum_cur[last])
    med = float(cum_base.quantile(0.5, axis=1).iloc[last])
    row = cum_base.iloc[last].dropna()
    return dict(wide=wide, observed=observed, base_years=base_years, base=base,
                cur=cur, days=days, cum_base=cum_base, cum_cur=cum_cur,
                upto=int(base.index.to_numpy()[last]), now=now, med=med,
                rank=int((row < now).sum()) + 1, n_years=len(row) + 1,
                pct=100 * (now - med) / med if med else np.nan)


def season_integral(sub: pd.DataFrame, verbose: bool = False) -> dict | None:
    """Season NDMI integral of the current year against the baseline median.

    Args:
        sub: Cleaned rows (see `prepare`) for ONE site and regime.
        verbose: Print why years or the whole series were dropped.

    Returns:
        dict with integral_pct, rank (1 = lowest), n_years, upto (last DOY
        integrated), n_days (days shared by all years), window, and
        n_baseline_years; None if the data cannot support the comparison.
    """
    p = integral_parts(sub, verbose=verbose)
    if p is None:
        return None
    return {"integral_pct": float(p["pct"]), "rank": p["rank"], "n_years": p["n_years"],
            "upto": p["upto"], "n_days": int(p["days"].sum()),
            "window": tuple(INTEGRATION_WINDOW),
            "n_baseline_years": len(p["base_years"])}


