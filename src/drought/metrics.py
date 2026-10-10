"""Analysis metrics.

Which metric to use is not a free choice here, and the reasoning is worth
keeping attached to the code:

  cut COUNT       Not comparable between years. Detection sensitivity rises
                  steeply with the number of clear observations (~1.3 cuts at
                  <15 images, ~2.9 at 35-45), and drought summers are sunny, so
                  2026 has roughly twice the images of 2024. Use `equalised`,
                  which thins every parcel-year to the same count.
  late SHARE      A ratio, so sparse sampling loses cuts from both halves of the
                  season roughly equally. Better, but still needs care.
  LAST CUT        The most robust: the detector needs to find ONE cut per
                  parcel-year, not all of them, and the observational bias runs
                  the right way — 2026 has the most images, so a late cut is
                  easier to detect there, not harder.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .config import CFG
from .detect import detect_all, detect_parcel_year, season_end


def _truncate(df, cfg):
    if cfg.detect.truncate_to_current:
        return df[df.doy <= season_end(df, cfg)]
    return df


def _subset(df, regime, cfg):
    lu = cfg.landuse
    codes = lu.intensive if regime == "intensive" else lu.extensive
    return df[df.lnf_code.isin(codes)] if "lnf_code" in df else df


def summarise(ev: pd.DataFrame) -> pd.DataFrame:
    """Per-year cut counts, with the observation diagnostics beside them.

    n_obs and max_gap belong in the same table as n_cuts: if one year had fewer
    observations, a lower count could be an observing artefact. Reporting them
    together is what makes the count interpretable.
    """
    return ev.groupby("year").agg(
        parcels=("n_cuts", "size"), mean_cuts=("n_cuts", "mean"),
        median_cuts=("n_cuts", "median"), first_cut=("first_cut", "median"),
        mean_obs=("n_obs", "mean"), max_gap=("max_gap", "median"),
        amp=("amp", "mean")).round(2)


def late_fraction(df: pd.DataFrame, regime="intensive", cfg=None,
                  split=None) -> tuple:
    """Share of detected cuts after `split`, overall and by observation density.

    The density table is printed so a reader can confirm the ratio is stable
    WITHIN a year across densities, which is what licenses comparing it between
    years.
    """
    C = cfg or CFG
    split = split or C.analysis.late_split_doy
    d = _truncate(_subset(df, regime, C), C)
    rows = []
    for (pid, yr), g in d.groupby(["parcel_id", "year"]):
        r = detect_parcel_year(g, C)
        cd = r["cut_dates"] or []
        if not cd or np.isnan(r["n_cuts"]):
            continue
        rows.append({"parcel_id": pid, "year": yr, "n_obs": r["n_obs"],
                     "frac_late": sum(1 for x in cd if x >= split) / len(cd)})
    t = pd.DataFrame(rows)
    overall = t.groupby("year").agg(parcel_years=("frac_late", "size"),
                                    frac_late=("frac_late", "mean"),
                                    mean_obs=("n_obs", "mean")).round(3)
    t["bin"] = pd.cut(t.n_obs, [0, 15, 25, 35, 50])
    by_density = t.pivot_table(index="bin", columns="year", values="frac_late",
                               aggfunc="mean", observed=True).round(3)
    counts = t.pivot_table(index="bin", columns="year", values="frac_late",
                           aggfunc="size", observed=True).fillna(0).astype(int)
    return overall, by_density, counts


def equalised(df: pd.DataFrame, regime="intensive", cfg=None, n_rep=10,
              seed=0, target=None) -> pd.DataFrame:
    """Cut counts with observation density equalised across years.

    Randomly thin every parcel-year to the SAME number of observations, detect,
    and repeat. The mean over repeats is the count each year would have earned
    under identical observing conditions.
    """
    C = cfg or CFG
    rng = np.random.default_rng(seed)
    d = _truncate(_subset(df, regime, C), C)
    d = d[d.doy >= C.detect.season[0]]
    per = d.groupby(["parcel_id", "year"]).size()
    target = target or int(per.groupby(level=1).median().min())
    print(f"  thinning every parcel-year to {target} observations, "
          f"{n_rep} repeats")
    out = []
    for _ in range(n_rep):
        keep = [g if len(g) <= target
                else g.sample(target, random_state=int(rng.integers(1 << 31)))
                for _, g in d.groupby(["parcel_id", "year"])]
        ev = detect_all(pd.concat(keep), C, verbose=False)
        out.append(ev.groupby("year").n_cuts.mean())
    res = pd.concat(out, axis=1)
    return pd.DataFrame({"mean_cuts": res.mean(axis=1).round(2),
                         "sd_over_repeats": res.std(axis=1).round(3)})


def last_cut(df: pd.DataFrame, regime="intensive", cfg=None) -> pd.DataFrame:
    """Date of the final detected cut, per parcel-year."""
    C = cfg or CFG
    d = _truncate(_subset(df, regime, C), C)
    rows = []
    for (pid, yr), g in d.groupby(["parcel_id", "year"]):
        r = detect_parcel_year(g, C)
        cd = r["cut_dates"] or []
        if not cd:
            continue
        rows.append({"parcel_id": pid, "year": yr, "last_cut": max(cd),
                     "first_cut": min(cd), "n_cuts": len(cd),
                     "n_obs": r["n_obs"]})
    return pd.DataFrame(rows)


def last_cut_paired(lc: pd.DataFrame, cfg=None, verbose=True) -> pd.DataFrame:
    """Each parcel against its OWN baseline mean.

    Pairing removes site, farm and management style in one step. A between-year
    mean over different parcels cannot, because the parcels observed in one year
    are not the parcels observed in another.
    """
    C = cfg or CFG
    year = C.project.current_year
    piv = lc.pivot_table(index="parcel_id", columns="year", values="last_cut")
    base = [c for c in piv.columns if c != year]
    out = pd.DataFrame({"baseline": piv[base].mean(axis=1),
                        "current": piv.get(year)}).dropna()
    out["shift_days"] = out.current - out.baseline
    if verbose:
        print(f"{len(out)} parcels observed in {year} and >=1 baseline year")
        print(f"  baseline mean last cut  DOY {out.baseline.mean():.1f}")
        print(f"  {year} mean last cut     DOY {out.current.mean():.1f}")
        print(f"  median shift            {out['shift_days'].median():+.1f} days")
        print(f"  parcels cutting earlier {100 * (out['shift_days'] < 0).mean():.0f}%")
    return out


def labels_paired(labels: dict, cfg=None, end=None) -> pd.DataFrame:
    """The same paired comparison, from HAND LABELS — no detector involved.

    This is the most defensible version of the result: it depends on no
    algorithm, only on reading the curves. Both years are truncated to the same
    day so the comparison is like for like.
    """
    C = cfg or CFG
    lo, hi = C.detect.season
    end = end or hi
    rows = []
    for key, v in labels.items():
        pid, yr = key.split("|")
        v = [d for d in v if lo <= d <= end]
        if v:
            rows.append({"parcel_id": pid, "year": int(yr), "n": len(v),
                         "first": min(v), "last": max(v)})
    t = pd.DataFrame(rows)
    piv = t.pivot_table(index="parcel_id", columns="year",
                        values=["n", "first", "last"]).dropna()
    return piv
