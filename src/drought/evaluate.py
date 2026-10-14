"""Scoring the detector against hand labels.

What is measured: whether the detector puts a cut where a human put one. Each
detection is matched to an unmatched label within a tolerance in days, so two
detections near one real cut give one hit and one false positive, not two hits.

What is NOT measured: truth. The labels are read off the same NDMI curves the
detector uses, so a cut that left no trace in NDMI is invisible to both. This is
agreement with an expert reading of the signal, which is the right target, but
it is not a farm diary.

And the metric that matters for the headline is not F1. The claim is a
comparison BETWEEN years, so what threatens it is error rates that DIFFER by
year — see `score_by_year`.
"""

from __future__ import annotations

import inspect
import itertools

import numpy as np
import pandas as pd

from .config import CFG
from .detect import detect_parcel_year


def _match(pred, truth, tol):
    tp = fp = 0
    um = list(truth)
    for p in pred:
        hit = [t for t in um if abs(t - p) <= tol]
        if hit:
            um.remove(min(hit, key=lambda t: abs(t - p))); tp += 1
        else:
            fp += 1
    return tp, fp, len(um)


def score(labels: dict, df: pd.DataFrame, cfg=None, detector=None,
          tol=None, verbose=True, **over) -> dict:
    C = cfg or CFG
    tol = tol or C.evaluate.match_tolerance_days
    detector = detector or detect_parcel_year
    tp = fp = fn = 0
    for key, truth in labels.items():
        pid, yr = key.split("|")
        g = df[(df.parcel_id == pid) & (df.year == int(yr))]
        if g.empty:
            continue
        pred = detector(g, C, **over)["cut_dates"] or []
        a, b, c = _match(pred, truth, tol)
        tp += a; fp += b; fn += c
    prec = tp / (tp + fp) if tp + fp else np.nan
    rec = tp / (tp + fn) if tp + fn else np.nan
    f1 = 2 * prec * rec / (prec + rec) if prec and rec else np.nan
    if verbose:
        print(f"  labelled cuts {tp + fn} | detected {tp + fp}")
        print(f"  precision {prec:.3f}  recall {rec:.3f}  F1 {f1:.3f}")
    return {"precision": round(prec, 3), "recall": round(rec, 3),
            "f1": round(f1, 3), "tp": tp, "fp": fp, "fn": fn}


def score_by_year(labels: dict, df: pd.DataFrame, cfg=None, tol=None,
                  detector=None) -> pd.DataFrame:
    """Precision and recall PER YEAR.

    Flat errors across years are what license a between-year comparison. A
    detector at F1 0.89 with year-neutral errors is more useful here than one
    at 0.94 whose false positives cluster in the drought year.
    """
    C = cfg or CFG
    tol = tol or C.evaluate.match_tolerance_days
    detector = detector or detect_parcel_year
    rows = []
    for key, truth in labels.items():
        pid, yr = key.split("|")
        g = df[(df.parcel_id == pid) & (df.year == int(yr))]
        if g.empty:
            continue
        pred = detector(g, C)["cut_dates"] or []
        tp, fp, fn = _match(pred, truth, tol)
        rows.append({"year": int(yr), "tp": tp, "fp": fp, "fn": fn})
    r = pd.DataFrame(rows)
    t = r.groupby("year").sum()
    t["precision"] = (t.tp / (t.tp + t.fp)).round(3)
    t["recall"] = (t.tp / (t.tp + t.fn)).round(3)
    t["n_parcel_years"] = r.groupby("year").size()
    return t


def miss_positions(labels: dict, df: pd.DataFrame, cfg=None, tol=None,
                   split=None) -> pd.DataFrame:
    """WHERE the misses fall in the season.

    Only misses that are a parcel's LAST cut can bias the last-cut metric. If
    the extra misses in one year are early-season, that metric is untouched.
    """
    C = cfg or CFG
    tol = tol or C.evaluate.match_tolerance_days
    split = split or C.analysis.late_split_doy
    rows = []
    for key, truth in labels.items():
        pid, yr = key.split("|")
        g = df[(df.parcel_id == pid) & (df.year == int(yr))]
        if g.empty:
            continue
        pred = detect_parcel_year(g, C)["cut_dates"] or []
        um = list(truth)
        for p in pred:
            hit = [t for t in um if abs(t - p) <= tol]
            if hit:
                um.remove(min(hit, key=lambda t: abs(t - p)))
        for m in um:
            rows.append({"year": int(yr), "doy": m, "key": key,
                         "is_last": m == max(truth),
                         "half": "late" if m >= split else "early"})
    return pd.DataFrame(rows)


def false_positive_positions(labels: dict, df: pd.DataFrame, cfg=None, tol=None
                             ) -> pd.DataFrame:
    """WHERE the false positives fall: one row per detection with no label in reach.

    A drought decline or its recovery can look like a cut. If the drought year's
    false positives sit in high summer and in the recovery, that is the
    mechanism; if they are spread like everyone else's, it is not.
    """
    C = cfg or CFG
    tol = tol or C.evaluate.match_tolerance_days
    rows = []
    for key, truth in labels.items():
        pid, yr = key.split("|")
        g = df[(df.parcel_id == pid) & (df.year == int(yr))]
        if g.empty:
            continue
        um = list(truth)
        for p in sorted(detect_parcel_year(g, C)["cut_dates"] or []):
            hit = [t for t in um if abs(t - p) <= tol]
            if hit:
                um.remove(min(hit, key=lambda t: abs(t - p)))
            else:
                rows.append({"year": int(yr), "doy": int(p), "key": key})
    return pd.DataFrame(rows, columns=["year", "doy", "key"])


def assert_wired(cfg=None) -> bool:
    """Check every detection parameter is actually READ by the detector.

    A parameter can exist in the config while the code that uses it is an older
    copy. Overriding it then changes nothing, and an A/B test reports "this rule
    does not help" when the truth is "this rule was never connected". That
    happened once and cost a wrong conclusion.
    """
    C = cfg or CFG
    src = inspect.getsource(detect_parcel_year)
    dead = [k for k in vars(C.detect)
            if k not in ("index", "confirm", "season", "truncate_to_current")
            and f'"{k}"' not in src and f"'{k}'" not in src]
    if dead:
        raise RuntimeError(
            f"these config parameters are NOT used by detect_parcel_year: "
            f"{dead}. The detector is out of date, or the names have drifted "
            "from the YAML — any comparison involving them is meaningless.")
    return True


def compare(labels: dict, df: pd.DataFrame, configs: dict, cfg=None,
            tol=None) -> pd.DataFrame:
    """Score several named parameter sets against the SAME labels."""
    C = cfg or CFG
    assert_wired(C)
    rows = []
    for name, over in configs.items():
        r = score(labels, df, C, tol=tol, verbose=False, **over)
        rows.append({"config": name, "detected": r["tp"] + r["fp"],
                     "hits": r["tp"], "false_pos": r["fp"], "missed": r["fn"],
                     "precision": r["precision"], "recall": r["recall"],
                     "f1": r["f1"], **over})
    return pd.DataFrame(rows).set_index("config")


def tune(labels: dict, df: pd.DataFrame, grid: dict | None = None, cfg=None,
         tol=None) -> pd.DataFrame:
    """Grid search against hand labels.

    Prefer a setting on a PLATEAU of good scores over the single best cell: a
    lone peak among two dozen combinations is usually fitted to labelling noise.
    """
    C = cfg or CFG
    grid = grid or {"min_drop_abs": [0.16, 0.20, 0.24],
                    "min_concentration": [0.0, 0.6],
                    "require_followup": [0, 20],
                    "lookback_days": [14, 26]}
    keys = list(grid)
    rows = []
    for combo in itertools.product(*(grid[k] for k in keys)):
        over = dict(zip(keys, combo))
        r = score(labels, df, C, tol=tol, verbose=False, **over)
        rows.append({**over, "precision": r["precision"],
                     "recall": r["recall"], "f1": r["f1"]})
    return pd.DataFrame(rows).sort_values("f1", ascending=False)


def late_miss_report(labels: dict, df: pd.DataFrame, cfg=None, tol=None,
                     split=None, n_show: int = 6):
    """Which parcel-years lost a LATE cut, and how far it was from anything found.

    A late miss is the only kind that can bias the last-cut metric, and it is
    also the kind that is easiest to fix once you can see the curve. This lists
    them with a `gap_days` column: the distance from the missed cut to the
    nearest detection. A large gap means the detector saw nothing there at all;
    a small one means it fired slightly off and the tolerance rejected the
    match, which is a different problem with a different fix.
    """
    C = cfg or CFG
    tol = tol or C.evaluate.match_tolerance_days
    split = split or C.analysis.late_split_doy
    rows = []
    for key, truth in labels.items():
        pid, yr = key.split("|")
        g = df[(df.parcel_id == pid) & (df.year == int(yr))]
        if g.empty:
            continue
        pred = detect_parcel_year(g, C)["cut_dates"] or []
        um = list(truth)
        for p in pred:
            hit = [t for t in um if abs(t - p) <= tol]
            if hit:
                um.remove(min(hit, key=lambda t: abs(t - p)))
        for m in um:
            if m < split:
                continue
            gap = min((abs(m - p) for p in pred), default=np.nan)
            n_obs = int(((g.doy >= m - 20) & (g.doy <= m + 20) &
                         (~g.haze if "haze" in g else True)).sum())
            rows.append({"key": key, "year": int(yr), "missed_doy": m,
                         "gap_days": gap, "obs_within_20d": n_obs,
                         "is_last": m == max(truth)})
    t = pd.DataFrame(rows)
    if t.empty:
        print("no late misses")
        return t
    print(f"{len(t)} late misses (DOY >= {split}), "
          f"{int(t.is_last.sum())} of them a parcel's LAST cut\n")
    print(t.groupby("year").agg(late_misses=("missed_doy", "size"),
                                were_last=("is_last", "sum"),
                                median_gap=("gap_days", "median"),
                                median_obs_nearby=("obs_within_20d", "median")
                                ).to_string())
    print("\n  gap_days: distance to the nearest detection. Large = nothing "
          "found there;\n  small = fired slightly off and the tolerance "
          "rejected it.")
    print("\nworst cases (inspect these with figures.plot_parcel):")
    print(t.nlargest(min(n_show, len(t)), "gap_days")
           [["key", "missed_doy", "gap_days", "obs_within_20d", "is_last"]]
           .to_string(index=False))
    return t


def score_by_half(labels: dict, df: pd.DataFrame, cfg=None, tol=None,
                  split=None, by_year: bool = False, **over) -> pd.DataFrame:
    """Precision and recall split into EARLY and LATE cuts.

    Why this and not overall F1: early cuts outnumber late ones, so the overall
    score is dominated by them. A detector can sit at F1 0.91 and still be poor
    at exactly the cuts the last-cut metric depends on. Tune against this, not
    against the aggregate.

    A detection is counted as early or late by the LABEL it matched, so a
    detection that matched nothing is attributed by its own date.
    """
    C = cfg or CFG
    tol = tol or C.evaluate.match_tolerance_days
    split = split or C.analysis.late_split_doy
    rows = []
    for key, truth in labels.items():
        pid, yr = key.split("|")
        g = df[(df.parcel_id == pid) & (df.year == int(yr))]
        if g.empty:
            continue
        pred = detect_parcel_year(g, C, **over)["cut_dates"] or []
        um = list(truth)
        for p in pred:
            hit = [t for t in um if abs(t - p) <= tol]
            if hit:
                m = min(hit, key=lambda t: abs(t - p))
                um.remove(m)
                rows.append({"year": int(yr), "half": "late" if m >= split
                             else "early", "kind": "tp"})
            else:
                rows.append({"year": int(yr), "half": "late" if p >= split
                             else "early", "kind": "fp"})
        for m in um:
            rows.append({"year": int(yr), "half": "late" if m >= split
                         else "early", "kind": "fn"})
    r = pd.DataFrame(rows)
    keys = ["half", "year"] if by_year else ["half"]
    t = (r.groupby(keys + ["kind"]).size().unstack("kind", fill_value=0)
          .reindex(columns=["tp", "fp", "fn"], fill_value=0))
    t["precision"] = (t.tp / (t.tp + t.fp)).round(3)
    t["recall"] = (t.tp / (t.tp + t.fn)).round(3)
    t["f1"] = (2 * t.precision * t.recall / (t.precision + t.recall)).round(3)
    return t


def tune_for(labels: dict, df: pd.DataFrame, grid: dict, cfg=None, tol=None,
             objective: str = "late_recall") -> pd.DataFrame:
    """Grid search against a chosen objective, not automatically overall F1.

    `objective` is one of: late_recall, late_f1, overall_f1, precision.

    Use late_recall when the question is missed late cuts — but read the
    precision column too. A setting that lifts late recall by inventing late
    cuts is worse than the problem, because a spurious late cut moves the
    last-cut metric in the direction you are trying to measure.
    """
    import itertools
    C = cfg or CFG
    assert_wired(C)
    keys = list(grid)
    rows = []
    for combo in itertools.product(*(grid[k] for k in keys)):
        over = dict(zip(keys, combo))
        h = score_by_half(labels, df, C, tol=tol, **over)
        o = score(labels, df, C, tol=tol, verbose=False, **over)
        row = {**over,
               "late_recall": h.loc["late", "recall"] if "late" in h.index else np.nan,
               "late_prec": h.loc["late", "precision"] if "late" in h.index else np.nan,
               "late_f1": h.loc["late", "f1"] if "late" in h.index else np.nan,
               "early_recall": h.loc["early", "recall"] if "early" in h.index else np.nan,
               "overall_f1": o["f1"], "precision": o["precision"]}
        rows.append(row)
    key = {"late_recall": "late_recall", "late_f1": "late_f1",
           "overall_f1": "overall_f1", "precision": "precision"}[objective]
    return pd.DataFrame(rows).sort_values(key, ascending=False)
