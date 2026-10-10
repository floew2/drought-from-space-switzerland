"""Synthetic Sentinel-2 index series: site trajectories and per-parcel cut series.

Two products, shaped exactly like the Earth Engine exports they stand in for:

trajectory   one row per site / regime / year / acquisition date, columns
             site, regime, year, date, doy, NDVI, kNDVI, NDMI, NDRE, n_px —
             as written by `gee/gee_grassland_trajectory.js`
cuts         one wide CSV per site and year, one row per parcel, one column per
             date and band (`<sensing>_<processing>_<tile>_<index>_<stat>`) — as
             written by `gee/gee_cuts_extract.js`, minus the geometry column

What is modelled, because the pipeline has to cope with it:

* cloud: dates are missing, and some are only partly clear; the drought summer
  is sunnier, so 2026 has more clear dates than the baseline years
* residual haze: NDVI pulled down while NDMI barely moves
* single-date dips that bounce straight back (shadow, misregistration)
* mowing: a step drop of 0.25-0.40 NDMI and regrowth over about three weeks;
  intensive meadows cut 4-6 times from early May, extensive ones about twice and
  not before 15 June
* drought: the 2026 canopy declines slowly through July (a slope, not a step),
  regrowth after a July cut stalls, late cuts are skipped, and the canopy
  recovers from early September

The site trajectories are CALIBRATED: the depth and onset of the 2026 collapse
are solved so that the trajectory module's own season integral reproduces the
real result (targets.INTEGRAL_PCT) for each site and regime.
"""

from __future__ import annotations

import datetime as dt

import numpy as np
import pandas as pd
from scipy import optimize

from . import targets as T

YEARS = list(range(2019, 2027))
CURRENT = 2026
DOY_RANGE = (60, 305)
END_2026 = {"linthebene": 248, "rheintal": 246, "alpstein": 250}

# Full-coverage pixel count of each site / regime mask (n_px on a clear date).
FULL_PX = {("linthebene", "intensive"): 6200, ("linthebene", "extensive"): 610,
           ("rheintal", "intensive"): 8100, ("rheintal", "extensive"): 720,
           ("alpstein", "intensive"): 240, ("alpstein", "extensive"): 380}
# Probability that a revisit is usable. Alpine cloud is the reason Alpstein
# keeps too few gap-free baseline years for a firm result.
P_CLEAR = {"linthebene": 0.62, "rheintal": 0.62, "alpstein": 0.36}
# Drought summers are sunny: more clear dates, so the deficit is measured under
# BETTER observing conditions than the baseline, not worse.
P_CLEAR_SUMMER = {2022: 0.75, 2026: 0.80}
TILES = ("T32TMT", "T32TNT")


def ndvi_from_ndmi(ndmi):
    """NDVI saturates over dense swards while NDMI keeps responding."""
    return 0.90 - 0.45 * np.exp(-5.0 * np.asarray(ndmi, float))


def ndre_from_ndmi(ndmi):
    return 0.12 + 0.85 * np.asarray(ndmi, float)


def _smoothstep(x, a, b):
    t = np.clip((np.asarray(x, float) - a) / (b - a), 0.0, 1.0)
    return t * t * (3 - 2 * t)


# --------------------------------------------------------------------------- #
# Acquisitions
# --------------------------------------------------------------------------- #
def acquisitions(rng: np.random.Generator, site: str, year: int) -> pd.DataFrame:
    """Usable acquisitions for one site and year: doy, coverage, haze, tile.

    Sentinel-2A/B revisit every 5 days, with a second orbit from the
    neighbouring track; a date is usable with a site- and season-dependent
    probability, and a quarter of usable dates are only partly clear.
    """
    phase = int(rng.integers(0, 5))
    doys = np.arange(DOY_RANGE[0] + phase, DOY_RANGE[1] + 1, 5)
    doys = np.sort(np.concatenate([doys, doys[rng.uniform(size=len(doys)) < 0.35] + 2]))
    p = np.full(len(doys), P_CLEAR[site])
    if year in P_CLEAR_SUMMER and site != "alpstein":
        p = np.where((doys >= 160) & (doys <= 245), P_CLEAR_SUMMER[year], p)
    if year == CURRENT:
        keep = doys <= END_2026[site]
        doys, p = doys[keep], p[keep]
        # the export date: the last acquisition before it is clear
        doys = np.append(doys, END_2026[site]) if doys[-1] != END_2026[site] else doys
        p = np.append(p, 1.0)[: len(doys)]
        p[-1] = 1.0
    use = rng.uniform(size=len(doys)) < p
    d = doys[use]
    cov = np.where(rng.uniform(size=len(d)) < 0.15, rng.uniform(0.2, 0.8, len(d)),
                   rng.uniform(0.9, 1.0, len(d)))
    if year == CURRENT:
        cov[-1] = max(cov[-1], 0.95)
    haze = rng.uniform(size=len(d)) < 0.07
    dip = (~haze) & (rng.uniform(size=len(d)) < 0.02)
    if year == CURRENT:
        haze[-1] = dip[-1] = False
    return pd.DataFrame({"doy": d, "cov": cov, "haze": haze, "dip": dip,
                         "tile": rng.choice(TILES, len(d))})


def _date(year: int, doy: int) -> dt.date:
    return dt.date(year, 1, 1) + dt.timedelta(days=int(doy) - 1)


# --------------------------------------------------------------------------- #
# Site trajectories
# --------------------------------------------------------------------------- #
def _baseline_ndmi(doy, regime, year_rng):
    """Site-mean NDMI for an ordinary year: green-up, mowing waves, autumn."""
    doy = np.asarray(doy, float)
    plateau = (0.36 if regime == "intensive" else 0.31) * year_rng["scale"]
    curve = 0.08 + (plateau - 0.08) * _smoothstep(doy, 75, 135)
    # Cuts are staggered across farms, so the site mean shows shallow waves,
    # not steps: the first-cut wave is the deepest.
    waves = ((140, 0.045), (182, 0.03), (220, 0.025)) if regime == "intensive" \
        else ((178, 0.04), (238, 0.025))
    for c, a in waves:
        curve -= a * np.exp(-0.5 * ((doy - c - year_rng["shift"]) / 9.0) ** 2)
    curve -= 0.10 * _smoothstep(doy, 265, 305)
    curve += year_rng["wiggle_a"] * np.sin(2 * np.pi * (doy - year_rng["wiggle_p"]) / 90)
    return curve


def _drought_2026(doy, k, minimum_drop):
    """Anomaly of 2026: slightly above normal in spring, collapse, recovery.

    The deficit reaches its full depth by mid-August in every case; k in
    [0, 1.4] sets how early the collapse starts (k = 0: early July, k = 1.4:
    late May). It is the one calibrated parameter, so it moves the integral
    without moving the August minimum. Recovery is sharp from early September,
    which the integral (ending DOY 245-250) only just reaches.
    """
    doy = np.asarray(doy, float)
    spring = 0.015 * _smoothstep(doy, 90, 120) * (1 - _smoothstep(doy, 150, 175))
    onset = 185.0 - 40.0 * k
    fall = minimum_drop * np.maximum(_smoothstep(doy, onset, onset + 45),
                                     _smoothstep(doy, 182, 228))
    recover = 1 - _smoothstep(doy, 245, 260)
    return spring - fall * recover


def _drought_2022(doy):
    doy = np.asarray(doy, float)
    return -0.04 * _smoothstep(doy, 185, 215) * (1 - _smoothstep(doy, 232, 250))


def _trajectory_rows(rng, site, regime, k, minimum_drop, acq_by_year, years_rng):
    rows = []
    full = FULL_PX[(site, regime)]
    for year in YEARS:
        acq = acq_by_year[year]
        base = _baseline_ndmi(acq.doy, regime, years_rng[year])
        if year == CURRENT:
            base = base + _drought_2026(acq.doy, k, minimum_drop)
        elif year == 2022:
            base = base + _drought_2022(acq.doy)
        ndmi = base + years_rng[year]["noise"][: len(acq)]
        ndvi = ndvi_from_ndmi(ndmi) + rng.normal(0, 0.008, len(acq))
        ndre = ndre_from_ndmi(ndmi) + rng.normal(0, 0.006, len(acq))
        h = acq.haze.to_numpy()
        ndvi = np.where(h, ndvi - rng.uniform(0.10, 0.20, len(acq)), ndvi)
        ndre = np.where(h, ndre - rng.uniform(0.04, 0.09, len(acq)), ndre)
        ndmi = np.where(h, ndmi - rng.uniform(0.0, 0.03, len(acq)), ndmi)
        for j, r in enumerate(acq.itertuples()):
            rows.append({"site": site, "regime": regime, "year": year,
                         "date": _date(year, r.doy).isoformat(), "doy": int(r.doy),
                         "NDVI": round(float(ndvi[j]), 5),
                         "kNDVI": round(float(np.tanh(ndvi[j] ** 2)), 5),
                         "NDMI": round(float(ndmi[j]), 5),
                         "NDRE": round(float(ndre[j]), 5),
                         "n_px": int(round(full * r.cov))})
    return rows


def _draw_ok(site, acq_by_year, setup, tj) -> bool:
    """Is this cloud draw consistent with what the real data looked like?

    The integral uses only days observed in every year, and in the real data
    that was 75-100 of the 153 days in the window, with all eight years kept.
    The cloud draw decides both, so a draw is accepted only if pilot series
    land in that range for both regimes — and the real integral is reachable
    within the collapse model's range.
    """
    lo, hi = T.SHARED_DAYS
    for regime, (years_rng, sub_seed) in setup.items():
        drop = _minimum_drop(site, regime, years_rng)
        res = {}
        for k in (0.0, 1.4):
            df = pd.DataFrame(_trajectory_rows(np.random.default_rng(sub_seed), site, regime,
                                               k, drop, acq_by_year, years_rng))
            res[k] = tj.season_integral(tj.prepare(df, verbose=False))
            if res[k] is None or not (lo <= res[k]["n_days"] <= hi) \
                    or res[k]["n_years"] != 8:
                return False
        target = T.INTEGRAL_PCT.get((site, regime))
        if target is not None and not (res[1.4]["integral_pct"] + 2 <= target
                                       <= res[0.0]["integral_pct"] - 2):
            return False
    return True


def _minimum_drop(site, regime, years_rng):
    """Deficit depth that puts the mid-August 2026 minimum inside the real range."""
    lo, hi = T.NDMI_MIN_AUGUST
    target = {"linthebene": lo + 0.25 * (hi - lo), "rheintal": lo + 0.75 * (hi - lo),
              "alpstein": hi}[site]
    if regime == "extensive":
        target -= 0.01
    base = float(_baseline_ndmi(228, regime, years_rng[CURRENT]))
    return max(base - target, 0.0)


def trajectories(rng: np.random.Generator, cfg=None, verbose: bool = True) -> pd.DataFrame:
    """All site / regime / year trajectories, calibrated to the real integrals.

    For each site and regime with a real target, the onset of the 2026
    collapse is solved by root-finding so that `trajectory.season_integral` — the
    same function that produces the reported figure — returns the target.

    Args:
        rng: Seeded generator.
        cfg: Config (default: the loaded one).
        verbose: Print the calibration result per site and regime.

    Returns:
        DataFrame in the trajectory export schema.
    """
    from .. import trajectory as tj

    out = []
    for site in ("linthebene", "rheintal", "alpstein"):
        calibrated = any(s_ == site for s_, _ in T.INTEGRAL_PCT)
        for attempt in range(100):
            n_draws = attempt + 1
            acq_by_year = {y: acquisitions(rng, site, y) for y in YEARS}
            setup = {}
            for regime in ("intensive", "extensive"):
                years_rng = {y: {"scale": rng.normal(1.0, 0.035), "shift": rng.normal(0, 5),
                                 "wiggle_a": rng.normal(0, 0.006),
                                 "wiggle_p": rng.uniform(0, 90),
                                 "noise": rng.normal(0, 0.010, 200)} for y in YEARS}
                setup[regime] = (years_rng, int(rng.integers(1 << 31)))
            if not calibrated or _draw_ok(site, acq_by_year, setup, tj):
                break
        else:
            raise RuntimeError(f"{site}: no cloud draw matched the observed record")
        if verbose and calibrated:
            print(f"    {site}: cloud draw {n_draws} accepted "
                  f"({T.SHARED_DAYS[0]}-{T.SHARED_DAYS[1]} shared days, 8 years, as observed)")
        for regime in ("intensive", "extensive"):
            years_rng, sub_seed = setup[regime]
            drop = _minimum_drop(site, regime, years_rng)
            target = T.INTEGRAL_PCT.get((site, regime))

            def build(k, drop=drop, sub_seed=sub_seed, acq_by_year=acq_by_year,
                      years_rng=years_rng, site=site, regime=regime):
                r = np.random.default_rng(sub_seed)
                return pd.DataFrame(_trajectory_rows(r, site, regime, k, drop,
                                                     acq_by_year, years_rng))

            if target is None:
                out.append(build(0.5))
                continue

            def pct(k, build=build):
                return tj.season_integral(tj.prepare(build(k), verbose=False))["integral_pct"]

            lo, hi = 0.0, 1.4                 # larger k = earlier collapse = deeper
            p_lo, p_hi = pct(lo), pct(hi)
            if not (p_hi <= target <= p_lo):
                raise RuntimeError(
                    f"{site}/{regime}: target {target}% outside the reachable "
                    f"range [{p_hi:.1f}, {p_lo:.1f}]% — revisit the collapse model")
            k = optimize.brentq(lambda k_, pct=pct, target=target: pct(k_) - target,
                                lo, hi, xtol=1e-4)
            df = build(k)
            if verbose:
                res = tj.season_integral(tj.prepare(df, verbose=False))
                print(f"    {site:<11}{regime:<10} k {k:.3f}  "
                      f"integral {res['integral_pct']:+.1f}% (target {target:+.0f}%)  "
                      f"rank {res['rank']}/{res['n_years']}  to DOY {res['upto']}  "
                      f"{res['n_days']} shared days")
            out.append(df)
    return pd.concat(out, ignore_index=True)


# --------------------------------------------------------------------------- #
# Per-parcel cut series
# --------------------------------------------------------------------------- #
def _cut_schedule(rng, regime, year):
    """True cut days for one parcel-year."""
    if regime == "intensive":
        days = [rng.normal(128, 6)]
        while True:
            nxt = days[-1] + rng.normal(38, 5)
            if nxt > 280:
                break
            days.append(nxt)
    else:
        days = [max(rng.normal(180, 8), 167.0), rng.normal(242, 10)]
        if rng.uniform() < 0.25:
            days.append(days[-1] + rng.normal(35, 4))
    days = [d for d in days if d <= 285]
    if year == CURRENT:
        # The drought: after mid-July there is almost nothing left to cut.
        days = [d for d in days if d < 196 or rng.uniform() < 0.15]
    elif year == 2022:
        days = [d for d in days if d < 205 or rng.uniform() < 0.6]
    return sorted(float(round(d)) for d in days)


def _parcel_curve(rng, regime, year, cuts, days):
    """Daily NDMI for one parcel-year, cuts and drought included."""
    plateau = (0.42 if regime == "intensive" else 0.36) + rng.normal(0, 0.025)
    env = 0.06 + (plateau - 0.06) * _smoothstep(days, 75, 125)
    env -= 0.12 * _smoothstep(days, 268, 305)
    if year == CURRENT:
        # A slope, not a step: ~0.3 over two months. The steepest point of a
        # smoothstep is 1.5x its mean slope, so the ramp spans 60 days to stay
        # under the detector's 0.010 / day drop-rate floor (max ~0.008 / day).
        env -= (env - 0.11) * _smoothstep(days, 175, 235) * (1 - _smoothstep(days, 244, 262))
    elif year == 2022:
        env -= 0.06 * _smoothstep(days, 190, 215) * (1 - _smoothstep(days, 232, 250))
    curve = env.copy()
    for c in cuts:
        # A cut removes the canopy at once; the stubble stays short for about a
        # week while the cut grass dries and is carted off, then the sward
        # regrows over roughly three weeks.
        depth = rng.uniform(0.25, 0.40)
        lag = rng.uniform(4.0, 8.0)
        grow = rng.uniform(18.0, 26.0)
        if year == CURRENT and 185 <= c < 244:
            grow *= 3.0                    # regrowth stalls in the drought
        dtt = days - c
        removed = depth * np.clip(1.0 - (dtt - lag) / grow, 0.0, 1.0)
        removed[dtt < 0] = 0.0
        curve = np.minimum(curve, env - removed)
    return np.clip(curve, 0.02, None)


def cut_series(rng: np.random.Generator, parcels: pd.DataFrame, sites=("linthebene", "rheintal"),
               n_per_site: int = 80):
    """Per-parcel series, labels and truth for the cut detector.

    Args:
        rng: Seeded generator.
        parcels: The synthetic 2025 parcels (from `parcels.generate`).
        sites: Sites to extract, as the real extraction does.
        n_per_site: Parcels sampled per site.

    Returns:
        (wide, truth): `wide` maps (site, year) -> DataFrame in the GEE export
        schema; `truth` has one row per parcel-year with the true cut days and
        the observations a human would see.
    """
    from ..config import CFG

    lu = CFG.landuse
    wide, truth = {}, []
    for site in sites:
        pool = parcels[(parcels.site == site) & parcels.status.eq("stable")
                       & parcels.code.isin(list(lu.intensive) + list(lu.extensive))]
        intensive = pool[pool.code.isin(lu.intensive)]
        extensive = pool[pool.code.isin(lu.extensive)]
        n_e = min(len(extensive), int(0.3 * n_per_site))
        pick = pd.concat([intensive.sample(n_per_site - n_e, random_state=int(rng.integers(1 << 31))),
                          extensive.sample(n_e, random_state=int(rng.integers(1 << 31)))])
        for year in YEARS:
            acq = acquisitions(rng, site, year)
            days = np.arange(1, 366, dtype=float)
            cols = {}
            meta = {"parcel_id": [], "lnf_code": [], "area_ha": [], "elev_m": []}
            for p in pick.itertuples():
                regime = "intensive" if p.code in lu.intensive else "extensive"
                cuts = _cut_schedule(rng, regime, year)
                curve = _parcel_curve(rng, regime, year, cuts, days)
                full = max(int(p.area_ha * 100 * 0.84), 4)
                pid = f"SYN-{int(p.objectid)}"
                meta["parcel_id"].append(pid)
                meta["lnf_code"].append(int(p.code))
                meta["area_ha"].append(round(float(p.area_ha), 2))
                meta["elev_m"].append(int(p.elev_m))
                seen = []
                for a in acq.itertuples():
                    # On a partly clear date this parcel is either clear or not.
                    if a.cov < 0.9 and rng.uniform() > a.cov:
                        continue
                    ndmi = curve[int(a.doy) - 1] + rng.normal(0, 0.015)
                    ndvi = ndvi_from_ndmi(ndmi) + rng.normal(0, 0.012)
                    if a.haze:
                        ndvi -= rng.uniform(0.10, 0.20)
                        ndmi -= rng.uniform(0.0, 0.03)
                    if a.dip:
                        ndvi -= 0.20
                        ndmi -= 0.25
                    frac = 1.0 if rng.uniform() > 0.08 else rng.uniform(0.2, 0.55)
                    n_px = max(int(round(full * frac)), 1)
                    d = _date(year, a.doy).strftime("%Y%m%d")
                    stem = f"{d}T102021_{d}T102559_{a.tile}"
                    for name, val in ((f"{stem}_NDVI_median", ndvi), (f"{stem}_NDMI_median", ndmi),
                                      (f"{stem}_NDVI_count", n_px), (f"{stem}_NDMI_count", n_px)):
                        cols.setdefault(name, {})[pid] = round(float(val), 5)
                    seen.append((int(a.doy), bool(a.haze or a.dip), frac >= 0.6))
                truth.append({"site": site, "parcel_id": pid, "year": year,
                              "regime": regime, "cuts": cuts, "seen": seen})
            df = pd.DataFrame(meta)
            band = pd.DataFrame(cols).reindex(df.parcel_id).reset_index(drop=True)
            wide[(site, year)] = pd.concat([df, band[sorted(band.columns)]], axis=1)
    return wide, pd.DataFrame(truth)


def hand_labels(rng: np.random.Generator, truth: pd.DataFrame, n_by_site=None) -> dict:
    """Simulated hand labels: what a careful human marks on the curves.

    A labeller sees only the observations. A cut is marked at the first clear,
    haze-free observation within 12 days after it — the trough, as the real
    labelling instructions require. A cut with no such observation left no
    visible trace and goes unlabelled, exactly as it would by eye.

    Args:
        rng: Seeded generator.
        truth: From `cut_series`.
        n_by_site: Parcel-years to label per site; defaults to 40 + 25 = 65,
            the size of the real reference set.

    Returns:
        {site: {"<parcel_id>|<year>": [doy, ...]}}
    """
    from ..config import CFG

    n_by_site = n_by_site or {"linthebene": 40, "rheintal": 25}
    lo, hi = CFG.detect.season
    out = {}
    for site, n in n_by_site.items():
        t = truth[truth.site == site]
        # Spread across years, so per-year recall is measurable.
        per_year = max(n // t.year.nunique(), 1)
        pick = pd.concat([g.sample(min(per_year, len(g)),
                                   random_state=int(rng.integers(1 << 31)))
                          for _, g in t.groupby("year")])
        rest = t.drop(pick.index)
        pick = pd.concat([pick, rest.sample(n - len(pick), random_state=int(rng.integers(1 << 31)))])
        labels = {}
        for r in pick.itertuples():
            obs = [d for d, bad, clear in r.seen if clear and not bad]
            marks = []
            for c in r.cuts:
                after = [d for d in obs if c <= d <= c + 12]
                if after and lo <= after[0] <= hi:
                    marks.append(after[0])
            labels[f"{r.parcel_id}|{r.year}"] = sorted(set(marks))
        out[site] = labels
    return out
