"""Cut detection — the heuristic.

FROZEN configuration. Tuned on 42 hand-labelled parcel-years, then validated on
all 65: F1 0.911 -> 0.912. Re-tuning on the full label set would forfeit that
out-of-sample property, so don't, without a reason.

The rules and why each exists, in the order they were forced on us:

  index NDMI, not NDVI  Over dense grassland NDVI SATURATES near 0.90 and falls
                        only ~0.18 at a cut, while NDMI falls 0.25-0.40. NDMI
                        is also far less sensitive to residual haze, because
                        both its bands sit at long wavelengths.
  drop size             a cut is a real fall, not noise
  drop RATE             mowing is a step, drought senescence is a slope. Without
                        this the detector walks down a decline and calls each
                        stage a cut.
  confirm index         a real cut removes biomass, so NDVI must fall too
  persistence           cut grass stays short for weeks; an instant bounce-back
                        is atmosphere
  follow-up             an event with nothing after it cannot be verified, and
                        the haze filter is least reliable at series edges
  regrowth path         a cut inside a data gap has no visible falling limb,
                        but its recovery is still there
  green-up guard        ...and spring green-up is ALSO a low point followed by a
                        steep rise, so require the canopy to have been high
                        earlier. Without this every parcel gained a spurious
                        March cut.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .config import CFG


def regime_of(code, cfg=None) -> str:
    lu = (cfg or CFG).landuse
    try:
        c = int(code)
    except (TypeError, ValueError):
        return "unknown"
    if c in lu.intensive:
        return "intensive"
    if c in lu.extensive:
        return "extensive"
    return "unknown"


def season_end(df: pd.DataFrame, cfg=None) -> int:
    """Last day of year the current year actually reaches."""
    c = cfg or CFG
    cur = df[df.year == c.project.current_year]
    if cur.empty:
        return c.detect.season[1]
    return int(min(c.detect.season[1], cur.doy.max()))


def detect_parcel_year(g: pd.DataFrame, cfg=None, **over) -> dict:
    """Cut events for one parcel in one season.

    Works on the OBSERVATIONS, not an interpolated grid: interpolating across a
    three-week gap would invent the very shape we are trying to detect.

    An event is dated at the FIRST OBSERVATION THAT IS CLEARLY DOWN — the
    trough. Directly observable, so a human labelling by eye and the detector
    use the same rule, and both detection paths agree. It sits a median of ~3
    days after the actual cut, identically in every year.
    """
    C = cfg or CFG
    d = C.detect
    p = {k: over.get(k, getattr(d, k)) for k in vars(d)}
    p.update({k: v for k, v in over.items() if k not in p})

    regime = regime_of(g["lnf_code"].iloc[0], C) if "lnf_code" in g else "unknown"
    if "haze" in g.columns:
        g = g[~g["haze"]]
    lo, hi = p["season"]
    g = g[(g.doy >= lo) & (g.doy <= hi)].sort_values("doy")
    empty = {"n_cuts": np.nan, "first_cut": np.nan, "n_obs": len(g),
             "max_gap": np.nan, "amp": np.nan, "thr": np.nan,
             "cut_dates": [], "regime": regime}
    if len(g) < 6:
        return empty

    doy = g["doy"].to_numpy(float)
    sig = g[p["index"]].to_numpy(float)       # primary signal
    con = g[p["confirm"]].to_numpy(float)     # corroboration

    amp = float(np.percentile(sig, 90) - np.percentile(sig, 10))
    amp_c = float(np.percentile(con, 90) - np.percentile(con, 10))
    thr_season = max(p["min_drop_abs"], p["drop_fraction"] * amp)

    def _thr_at(i):
        """Detection threshold for the candidate at index i.

        With amp_window = 0 this is the season-wide threshold, unchanged. With a
        window, the amplitude is measured LOCALLY: late in a dry season the
        canopy is already low, so a real cut produces a smaller absolute drop
        than the same cut in May, and a threshold set by the whole-season
        amplitude is too high to catch it.
        """
        w = p.get("amp_window", 0)
        if not w:
            return thr_season
        near = (doy >= doy[i] - w) & (doy <= doy[i] + w)
        if near.sum() < 4:
            return thr_season
        loc = sig[near]
        a = float(np.percentile(loc, 90) - np.percentile(loc, 10))
        return max(p["min_drop_abs"], p["drop_fraction"] * a)

    events = []
    for i in range(1, len(doy)):
        back = (doy >= doy[i] - p["lookback_days"]) & (doy < doy[i])
        if not back.any():
            continue
        j = int(np.where(back)[0][np.argmax(sig[back])])
        drop = sig[j] - sig[i]
        thr = _thr_at(i)

        if drop < thr:
            # No usable falling limb — the cut may sit inside a gap. Look
            # forward: a steep recovery out of this trough implies biomass was
            # removed shortly before it.
            fwd = (doy > doy[i]) & (doy <= doy[i] + p["rise_window"])
            if not fwd.any():
                continue
            k = int(np.where(fwd)[0][np.argmax(sig[fwd])])
            rise = sig[k] - sig[i]
            if rise < thr or rise / max(doy[k] - doy[i], 1.0) < p["min_rise_rate"]:
                continue
            if sig[i] > np.percentile(sig, 35):
                continue
            prior = sig[:i]                   # green-up guard
            if prior.size == 0 or prior.max() - sig[i] < thr:
                continue
            events.append(doy[i])
            continue

        span = max(doy[i] - doy[j], 1.0)
        if drop / span < p["min_drop_rate"]:
            continue

        if p["min_concentration"] and i - j >= 2:
            steps = sig[j:i] - sig[j + 1:i + 1]
            if steps.max() < p["min_concentration"] * drop:
                continue

        dc = con[j] - con[i]
        if amp_c > 0 and (dc / amp_c) < p["confirm_frac"] * (drop / max(amp, 1e-6)):
            continue

        if i + 1 < len(doy) and (p["persist_always"]
                                 or (doy[i + 1] - doy[i]) <= 12):
            if sig[i + 1] > sig[j] - p["persist_frac"] * drop:
                continue

        if p["require_followup"]:
            if not ((doy > doy[i]) & (doy <= doy[i] + p["require_followup"])).any():
                continue

        events.append(doy[i])

    merged = []
    for e in sorted(events):
        if not merged or e - merged[-1] > p["merge_days"]:
            merged.append(e)

    return {"n_cuts": len(merged),
            "first_cut": merged[0] if merged else np.nan,
            "last_cut": merged[-1] if merged else np.nan,
            "n_obs": len(g),
            "max_gap": float(np.diff(doy).max()) if len(doy) > 1 else np.nan,
            "amp": amp, "thr": thr_season, "cut_dates": merged,
            "regime": regime}


def detect_all(df: pd.DataFrame, cfg=None, verbose: bool = True,
               **over) -> pd.DataFrame:
    C = cfg or CFG
    if C.detect.truncate_to_current:
        end = season_end(df, C)
        if end < C.detect.season[1] and verbose:
            print(f"  truncating ALL years at DOY {end} "
                  f"({C.project.current_year}'s last observation) so the "
                  "comparison is like for like")
        df = df[df.doy <= end]
    rows = []
    for (pid, yr), g in df.groupby(["parcel_id", "year"]):
        r = detect_parcel_year(g, C, **over)
        r.update(parcel_id=pid, year=yr,
                 elev_m=g["elev_m"].iloc[0] if "elev_m" in g else np.nan,
                 area_ha=g["area_ha"].iloc[0] if "area_ha" in g else np.nan,
                 lnf_code=g["lnf_code"].iloc[0] if "lnf_code" in g else np.nan)
        rows.append(r)
    return pd.DataFrame(rows)
