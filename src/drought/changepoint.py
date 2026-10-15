"""Changepoint detection — the alternative that was tested and rejected.

A penalised piecewise-linear fit: find the segmentation that best explains the
season, paying a penalty per break so it cannot simply split everywhere. Each
heuristic rule becomes a consequence rather than a rule — a slow drought
decline is fitted BY a sloping segment, so no break is placed there.

Scored on the same 65 hand-labelled parcel-years: best F1 0.871 against the
heuristic's 0.911, over-segmenting (32 false positives against 12); on the
10 m mask and 136 parcel-years, 0.805 against 0.868. NOT
ADOPTED. Kept as the record of a measured negative result, which is worth more
than an untested assumption in either direction.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .config import CFG
from .detect import detect_parcel_year, regime_of


def _segment_cost(t, y, a, b):
    """Residual sum of squares of a straight line fitted to points a..b."""
    n = b - a + 1
    if n <= 2:
        return 0.0
    tt, yy = t[a:b + 1], y[a:b + 1]
    tm, ym = tt.mean(), yy.mean()
    dt = tt - tm
    denom = float((dt * dt).sum())
    slope = float((dt * (yy - ym)).sum()) / denom if denom > 1e-9 else 0.0
    resid = yy - (ym + slope * dt)
    return float((resid * resid).sum())


def detect_changepoints(g, penalty_k=None, min_seg_days=None,
                        min_jump=None, index=None):
    """Cut dates from a penalised piecewise-linear fit.

    Returns the same shape as detect_parcel_year, so the two are interchangeable
    downstream and can be scored against the same labels.
    """
    penalty_k = CFG.changepoint.penalty_k if penalty_k is None else penalty_k
    min_seg_days = CFG.changepoint.min_seg_days if min_seg_days is None else min_seg_days
    min_jump = CFG.changepoint.min_jump if min_jump is None else min_jump
    index = CFG.detect.index if index is None else index

    regime = regime_of(g["lnf_code"].iloc[0]) if "lnf_code" in g else "unknown"
    if "haze" in g.columns:
        g = g[~g["haze"]]
    g = g[(g.doy >= CFG.detect.season[0]) & (g.doy <= CFG.detect.season[1])].sort_values("doy")
    empty = {"n_cuts": np.nan, "first_cut": np.nan, "n_obs": len(g),
             "max_gap": np.nan, "amp": np.nan, "thr": np.nan,
             "cut_dates": [], "regime": regime}
    if len(g) < 8:
        return empty

    t = g["doy"].to_numpy(float)
    y = g[index].to_numpy(float)
    n = len(t)

    # Noise scale from successive differences: robust to the sawtooth itself,
    # because a median absorbs the few large steps that are cuts.
    d = np.diff(y)
    sigma = float(np.median(np.abs(d - np.median(d))) * 1.4826 / np.sqrt(2))
    sigma = max(sigma, 0.005)
    penalty = penalty_k * sigma ** 2 * np.log(n)

    # A break is only CONSIDERED where the series falls. Without this the model
    # would happily place breaks at upward kinks, which are regrowth, not cuts.
    allowed = np.zeros(n, dtype=bool)
    for i in range(1, n):
        back = (t >= t[i] - CFG.detect.lookback_days) & (t < t[i])
        if back.any() and y[back].max() - y[i] >= min_jump:
            allowed[i] = True

    INF = float("inf")
    F = np.full(n + 1, INF)
    F[0] = -penalty          # the first segment is not charged a break
    prev = np.full(n + 1, -1, dtype=int)
    for j in range(1, n + 1):
        for i in range(0, j):
            if i > 0 and not allowed[i]:
                continue
            if t[j - 1] - t[i] < min_seg_days and j != n:
                continue
            if F[i] == INF:
                continue
            c = F[i] + _segment_cost(t, y, i, j - 1) + penalty
            if c < F[j]:
                F[j] = c
                prev[j] = i

    if F[n] == INF:
        return empty
    cuts, j = [], n
    while j > 0:
        i = prev[j]
        if i > 0:
            cuts.append(float(t[i]))
        j = i
    cuts.sort()

    return {"n_cuts": len(cuts),
            "first_cut": cuts[0] if cuts else np.nan,
            "n_obs": n,
            "max_gap": float(np.diff(t).max()),
            "amp": float(np.percentile(y, 90) - np.percentile(y, 10)),
            "thr": float(penalty),
            "cut_dates": cuts,
            "regime": regime}


def score_detector(labels, df, detector, tol=12, **kw):
    """Precision / recall for ANY detector with the standard return shape."""
    tp = fp = fn = 0
    for key, truth in labels.items():
        pid, yr = key.split("|")
        g = df[(df.parcel_id == pid) & (df.year == int(yr))]
        if g.empty:
            continue
        pred = detector(g, **kw)["cut_dates"] or []
        um = list(truth)
        for p in pred:
            hit = [x for x in um if abs(x - p) <= tol]
            if hit:
                um.remove(min(hit, key=lambda x: abs(x - p))); tp += 1
            else:
                fp += 1
        fn += len(um)
    prec = tp / (tp + fp) if tp + fp else np.nan
    rec = tp / (tp + fn) if tp + fn else np.nan
    return {"detected": tp + fp, "hits": tp, "false_pos": fp, "missed": fn,
            "precision": round(prec, 3), "recall": round(rec, 3),
            "f1": round(2 * prec * rec / (prec + rec), 3) if prec and rec
            else np.nan}


def tune_changepoint(labels, df, ks=(0.5, 1, 2, 3, 5, 8), tol=12):
    """Sweep the single sensitivity parameter, and show the heuristic beside it.

    One parameter to defend instead of six. Prefer a plateau: if only one value
    of k scores well, it is fitted to 114 labels rather than to the phenomenon.
    """
    rows = []
    for k in ks:
        r = score_detector(labels, df, detect_changepoints, tol=tol,
                           penalty_k=k)
        r["method"] = f"changepoint k={k}"
        rows.append(r)
    r = score_detector(labels, df, detect_parcel_year, tol=tol)
    r["method"] = "heuristic (6 rules)"
    rows.append(r)
    return pd.DataFrame(rows).set_index("method")[
        ["detected", "hits", "false_pos", "missed", "precision", "recall", "f1"]]


def plot_changepoint(df, pid, year, penalty_k=None):
    """One parcel-year with the fitted segments drawn, to see WHY it broke there."""
    g = df[(df.parcel_id == pid) & (df.year == year)]
    if "haze" in g:
        g = g[~g.haze]
    g = g[(g.doy >= CFG.detect.season[0]) & (g.doy <= CFG.detect.season[1])].sort_values("doy")
    r = detect_changepoints(g, penalty_k=penalty_k)
    t, y = g.doy.to_numpy(float), g[CFG.detect.index].to_numpy(float)

    plt.style.use("default")
    fig, ax = plt.subplots(figsize=(11, 3.2), dpi=130)
    ax.plot(t, y, "o", ms=4, color="#2166ac", label=CFG.detect.index)
    edges = [0] + [int(np.where(t == c)[0][0]) for c in r["cut_dates"]] + [len(t)]
    for a, b in zip(edges[:-1], edges[1:]):
        if b - a < 2:
            continue
        tt, yy = t[a:b], y[a:b]
        k = np.polyfit(tt, yy, 1)
        ax.plot(tt, np.polyval(k, tt), "-", color="#b2182b", lw=2)
    for c in r["cut_dates"]:
        ax.axvline(c, color="#ff3b30", lw=1.4, alpha=0.7)
    ax.set_title(f"{str(pid)[-8:]} {year} — {r['n_cuts']} breaks, "
                 f"penalty {r['thr']:.4f}", fontsize=9, loc="left")
    ax.set_xlabel("day of year"); ax.grid(alpha=0.3, ls=":")
    ax.legend(fontsize=8)
    plt.tight_layout(); plt.show()
    return r
