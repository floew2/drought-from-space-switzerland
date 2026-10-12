"""
Grassland VI trajectory and season integral — plots matching the precipitation
figure.

Input:  grassland_trajectory_NDVI.csv from gee_grassland_trajectory.js
        columns: site, regime, year, date, doy, ix_mean, n_px
Output: one two-panel figure per site/regime, plus a printed summary table.

Panel A: seasonal VI trajectory. Shows WHEN 2026 diverged, which a single
         anomaly number cannot.
Panel B: cumulative integral of VI. Direct analogue of the cumulative
         precipitation figure, and integrated VI is an established proxy for
         season biomass production.

Same visual grammar as precip_stgallen.py: percentile envelope, white dashed
median, red current year, square panels, dark theme.

Requires: pandas, numpy, matplotlib
"""

from __future__ import annotations

import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
from .config import CFG, figures_dir, gee_dir

CSV = None          # set by plot(); or pass a path directly

# Which index to plot. The CSV carries all four; switching is free.
#   kNDVI  recommended default — same information as NDVI but without the
#          saturation that compresses dense Swiss swards.
#   NDMI   canopy water; earlier signal, but snow inflates it (alpine caution).
#   NDRE   red-edge chlorophyll; Sentinel-2's real advantage over MODIS.
#   NDVI   familiar, weakest discriminator here.
INDEX = CFG.detect.index
ALL_INDICES = ["NDVI", "kNDVI", "NDMI", "NDRE"]
COMPARE_INDICES = True        # also print a side-by-side table of all four
SMOOTH_DAYS = CFG.analysis.smooth_days              # rolling-median window on the daily grid, 0 = off.
                              # A single hazy or partly-covered date should not
                              # be readable as a real drop; the raw observations
                              # stay visible as dots so nothing is hidden.
DIAGNOSE_COVERAGE = True      # report whether index value depends on coverage

# --- residual-haze filter --------------------------------------------------- #
# Cloud Score+ at cs_cdf >= 0.60 lets some thin haze and cirrus through. Its
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
HAZE_WINDOW = 5      # neighbouring observations used for the local median
HAZE_NDVI_Z = -1.5   # how far NDVI must fall below its local level (robust z)
HAZE_GAP = 1.5       # how much less NDMI fell, in the same units
CURRENT_YEAR = CFG.project.current_year
BASELINE = (CFG.project.baseline_years[0], CFG.project.baseline_years[-1])

DOY_MIN, DOY_MAX = 60, 305    # must match the GEE window

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

SITE_LABELS = {
    "rheintal": "St. Galler Rheintal",
    "linthebene": "Linthebene",
    "alpstein": "Alpstein / Churfirsten",
}
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
def load(path: str, verbose: bool = True) -> pd.DataFrame:
    """Read a trajectory export and clean it (see `prepare`)."""
    return prepare(pd.read_csv(path), verbose=verbose)


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
            print(f"  haze filter removed {n_haze} of {len(df)} dates "
                  f"({100 * n_haze / len(df):.1f}%)")
        if verbose and len(by_year):
            print("    by year: " + ", ".join(f"{int(y)}:{n}"
                                              for y, n in by_year.items())
                  + "  <- roughly even means the filter is not biasing a year")
        df = df[~df["haze"]].drop(columns="haze")

    df["ix_mean"] = df[INDEX]

    # Coverage relative to the best observation for that site and regime.
    full = df.groupby(["site", "regime"])["n_px"].transform("max")
    df["cov"] = df["n_px"] / full

    before = len(df)
    df = df[df["cov"] >= MIN_COVERAGE]
    if verbose:
        print(f"  dropped {before - len(df)} of {before} dates below "
              f"{MIN_COVERAGE:.0%} mask coverage")
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
                print(f"    [!] {int(year)}: only {len(g)} clear dates — skipped")
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


# --------------------------------------------------------------------------- #
# Plotting
# --------------------------------------------------------------------------- #
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


def _integral_parts(sub: pd.DataFrame, verbose: bool = True) -> dict | None:
    """Everything the season integral needs, shared by the number and the figure."""
    wide, observed = to_grid(sub, verbose=verbose)
    base_years = [y for y in wide.columns if BASELINE[0] <= y <= BASELINE[1]]
    if CURRENT_YEAR not in wide.columns or len(base_years) < 3:
        if verbose:
            print("    [!] insufficient data — skipped")
        return None

    lo, hi = INTEGRATION_WINDOW
    cov = window_coverage(wide[base_years + [CURRENT_YEAR]], lo, hi)
    # The current year is partial by definition — the season is not over. That
    # is not a data-quality problem and must not be reported as one; the
    # comparison simply ends where 2026 ends.
    thin = cov[base_years][cov[base_years] < MIN_WINDOW_COVERAGE]
    if verbose and not thin.empty:
        print("    [!] gappy baseline years dropped: "
              + ", ".join(f"{int(y)}:{c:.0%}" for y, c in thin.items()))
    base_years = [y for y in base_years if cov[y] >= MIN_WINDOW_COVERAGE]
    if len(base_years) < 3:
        if verbose:
            print("    [!] fewer than 3 usable baseline years — skipped")
        return None
    if verbose and len(base_years) < 5:
        print(f"    [!] only {len(base_years)} baseline years survive — the "
              "percentile envelope is weak here, treat the rank with caution")
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
    p = _integral_parts(sub, verbose=verbose)
    if p is None:
        return None
    return {"integral_pct": float(p["pct"]), "rank": p["rank"], "n_years": p["n_years"],
            "upto": p["upto"], "n_days": int(p["days"].sum()),
            "window": tuple(INTEGRATION_WINDOW),
            "n_baseline_years": len(p["base_years"])}


def make_figure(site: str, regime: str, sub: pd.DataFrame):
    p = _integral_parts(sub)
    if p is None:
        print(f"    [!] insufficient data for {site}/{regime} — skipped")
        return None
    lo, hi = INTEGRATION_WINDOW
    base, base_years, days = p["base"], p["base_years"], p["days"]
    cum_base, cum_cur, cur = p["cum_base"], p["cum_cur"], p["cur"]
    n_days_used = int(days.sum())
    print(f"    window DOY {lo}-{hi}, {len(base_years)} baseline years, "
          f"{CURRENT_YEAR} runs to DOY "
          f"{int(p['wide'].index[np.where(~np.isnan(cur))[0][-1]])}")
    print(f"    integral over {n_days_used} days shared by all "
          f"{len(base_years) + 1} years")
    if n_days_used < 40:
        print("    [!] few shared days — the integral is thin, prefer the "
              "trajectory panel and the mean anomaly")

    plt.style.use("dark_background")
    fig, axes = plt.subplots(2, 1, figsize=(6.4, 12.4), dpi=180,
                             constrained_layout=True)

    title = f"{SITE_LABELS.get(site, site)} - {REGIME_LABELS.get(regime, regime)}"
    obs_doy = sub.loc[sub["year"] == CURRENT_YEAR, "doy"].to_numpy()
    obs_val = sub.loc[sub["year"] == CURRENT_YEAR, INDEX].to_numpy()
    base_s = base.apply(lambda c: smooth(c.to_numpy(), SMOOTH_DAYS))
    base_s.index = base.index
    panel(axes[0], base_s, smooth(cur, SMOOTH_DAYS), INDEX, first=True,
          raw=(obs_doy, obs_val))
    axes[0].set_title("\n".join(textwrap.wrap(title, 34)), fontsize=12,
                      fontweight="bold", color="white", pad=10)
    panel(axes[1], cum_base, cum_cur,
          f"Cumulative {INDEX} (season integral)", first=False)
    axes[1].set_xlabel("Month", fontsize=10, color="white")

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

    fig.text(0.99, 0.004, f"Sentinel-2 / Cloud Score+ · {INDEX} · "
             "parcels: Nutzungsflächen",
             ha="right", va="bottom", fontsize=7, color="#6e7681")
    out = str(figures_dir(site) / f"vi_trajectory_{INDEX}_{regime}.png")
    fig.savefig(out, dpi=180, facecolor=fig.get_facecolor())
    print(f"  -> {out}")
    plt.close(fig)

    n_obs = {int(y): len(o) for y, o in p["observed"].items()}
    return dict(site=site, regime=regime, integral_pct=pct, rank=rank,
                upto=upto, window=(lo, hi), n_days=n_days_used,
                n_baseline_years=len(base_years), obs=n_obs, file=out)


# --------------------------------------------------------------------------- #
def _find_csv(site=None):
    """Locate the trajectory export.

    It covers ALL sites in one file, because the GEE script runs them together,
    so look in the site folder first (where you may have put it) and then in
    data/gee/ (where a multi-site file more naturally belongs).
    """
    from .config import ROOT
    cands = []
    if site:
        cands += sorted(gee_dir(site).glob("*trajectory*.csv"))
    cands += sorted((ROOT / CFG.paths.gee).glob("*trajectory*.csv"))
    if not cands:
        where = f"{gee_dir(site)} or {ROOT / CFG.paths.gee}" if site else str(ROOT / CFG.paths.gee)
        raise FileNotFoundError(f"no *trajectory*.csv in {where}")
    return cands[0]


def plot(site=None, path=None, index=None, regimes=None, verbose=True):
    """Trajectory and season-integral figures, one per site and regime.

    The trajectory panel is the trustworthy one: it shows WHEN the divergence
    happened and depends on no interpolation choices. The integral is sensitive
    to how gaps are handled, which is why it is the second panel and not the
    headline.
    """
    global INDEX
    if index:
        INDEX = index
    src = path or _find_csv(site)
    if verbose:
        print(f"reading {src}")
    from .provenance import verify
    verify(src, "trajectory")            # refuses a stale export, naming why
    df = load(str(src), verbose=verbose)   # this module's load() cleans inline

    sites = [site] if site else sorted(df.site.unique())
    out = []
    for s_ in sites:
        sub_all = df[df.site == s_]
        if sub_all.empty:
            print(f"  [!] no rows for site {s_!r} in this export")
            continue
        for regime, g in sub_all.groupby("regime"):
            if regimes and regime not in regimes:
                continue
            print(f"\n{s_} / {regime}: {len(g)} obs")
            r = make_figure(s_, regime, g)
            if r:
                out.append(r)
    return pd.DataFrame(out) if out else None


def main():
    df = load(CSV)
    print(f"{len(df):,} usable rows, sites: {sorted(df['site'].unique())}\n")

    results = []
    for (site, regime), sub in df.groupby(["site", "regime"]):
        print(f"{site} / {regime}: {len(sub)} obs")
        r = make_figure(site, regime, sub)
        if r:
            results.append(r)
            print(f"    integral to DOY {r['upto']}: {r['integral_pct']:+.1f}% "
                  f"vs. median, rank {r['rank']} -> {r['file']}")

    # Does the measured value depend on how much of the site was visible? If
    # low-coverage dates sit systematically lower, the jumps are an artefact of
    # partial scenes sampling different parcels, not real vegetation change.
    if DIAGNOSE_COVERAGE:
        print("\nCoverage diagnostic — mean index by coverage band:")
        print(f"  {'site':<12} {'regime':<7} {'0.8-0.9':>9} {'0.9-1.0':>9} {'r':>7}")
        for (site, regime), sub in df.groupby(["site", "regime"]):
            lo_b = sub[sub["cov"] < 0.9][INDEX].mean()
            hi_b = sub[sub["cov"] >= 0.9][INDEX].mean()
            r = sub[["cov", INDEX]].corr().iloc[0, 1]
            print(f"  {site:<12} {regime:<7} {lo_b:>9.3f} {hi_b:>9.3f} {r:>7.2f}")
        print("  a large gap between bands, or |r| > 0.3, means partial scenes "
              "are biasing the curve — raise MIN_COVERAGE.")

    # Index cross-comparison. If all four agree on sign and rough magnitude the
    # result is robust to index choice. Where they disagree, that is diagnostic:
    # NDMI diverging upward from NDVI at the alpine site points to snow, since
    # snow raises NDMI and lowers NDVI.
    if COMPARE_INDICES:
        present = [c for c in ALL_INDICES if c in df.columns]
        g = df[df["doy"].between(121, 243)]
        base_yrs = range(BASELINE[0], BASELINE[1] + 1)
        print("\nMay-Aug mean by index (2026 vs. baseline median of years):")
        hdr = "  ".join(f"{c:>8}" for c in present)
        print(f"  {'site':<12} {'regime':<7} {hdr}")
        for (site, regime), sub in g.groupby(["site", "regime"]):
            cur_ = sub[sub["year"] == CURRENT_YEAR]
            bas = sub[sub["year"].isin(base_yrs)]
            if cur_.empty or bas.empty:
                continue
            cells = []
            for c in present:
                b = bas.groupby("year")[c].mean().median()
                cells.append(f"{cur_[c].mean() - b:+8.3f}")
            print(f"  {site:<12} {regime:<7} " + "  ".join(cells))

    # Clear-observation diagnostic. Drought summers are sunny: if 2026 has MORE
    # clear dates than the baseline years, the deficit was measured under BETTER
    # observing conditions, which makes the finding conservative. Put this table
    # in the appendix — it converts an obvious objection into evidence.
    print("\nClear observations per year (site / regime):")
    for r in results:
        years = sorted(r["obs"])
        line = "  ".join(f"{y}:{r['obs'][y]}" for y in years)
        print(f"  {r['site']:<12} {r['regime']:<7} {line}")





def compare_masks(site=None, paths=None, index=None, split_doy=190):
    """Compare trajectories built from different mask edge treatments.

    The question: does tightening the purity threshold change the ANSWER, or
    only the noise? If a jump from 90% to 100% purity leaves the anomaly
    unmoved, the simple box purity is good enough and modelling the sensor
    point-spread function would be refining a smaller term than the LV95/UTM
    grid offset that GEE resolves on the fly anyway.

    `paths` maps a label to a trajectory CSV, e.g.
        {"purity 0.90": "data/gee/trajectory_multi_p090.csv",
         "purity 1.00": "data/gee/trajectory_multi_p100.csv"}
    """
    from .config import CFG as _C
    from .config import ROOT
    global INDEX
    if index:
        INDEX = index
    if paths is None:
        found = sorted((ROOT / _C.paths.gee).glob("*trajectory*p*.csv"))
        if len(found) < 2:
            raise FileNotFoundError(
                "need at least two trajectory exports built from different "
                f"mask bands; found {[f.name for f in found]}")
        paths = {f.stem.split("_")[-1]: f for f in found}

    cur = _C.project.current_year
    rows = []
    for label, path in paths.items():
        df = load(str(path))
        if site:
            df = df[df.site == site]
        g = df[df.doy.between(121, 243)]
        for (s_, reg), sub in g.groupby(["site", "regime"]):
            piv = sub.groupby("year")[INDEX].mean()
            base = piv[[y for y in piv.index if y != cur]]
            if cur not in piv.index or base.empty:
                continue
            rows.append({"mask": label, "site": s_, "regime": reg,
                         "n_obs": len(sub),
                         f"{cur}": round(float(piv[cur]), 4),
                         "baseline_median": round(float(base.median()), 4),
                         "anomaly": round(float(piv[cur] - base.median()), 4),
                         "rel_%": round(100 * float(piv[cur] - base.median())
                                        / float(base.median()), 1)})
    t = pd.DataFrame(rows)
    if t.empty:
        print("nothing to compare")
        return t
    print(f"May-Aug mean {INDEX}, {cur} against the baseline median:\n")
    print(t.to_string(index=False))
    piv = t.pivot_table(index=["site", "regime"], columns="mask",
                        values="anomaly")
    if piv.shape[1] > 1:
        spread = (piv.max(axis=1) - piv.min(axis=1)).abs()
        print(f"\nlargest difference between mask treatments: "
              f"{spread.max():.4f} {INDEX}")
        print(f"typical anomaly being measured:              "
              f"{t.anomaly.abs().median():.4f} {INDEX}")
        r = spread.max() / max(t.anomaly.abs().median(), 1e-9)
        print(f"\nratio {r:.1%} — under ~10% means the mask edge treatment is "
              "not\nwhat limits this result, and PSF modelling would refine a "
              "smaller term still.")
    return t
