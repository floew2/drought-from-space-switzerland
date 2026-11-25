"""Purity-threshold sensitivity figure: does the answer depend on the threshold?

The parcels that enter are those with at least one cell at
sensitivity.strict_purity. At the default 0.0 that is every parcel: the sweep
covers the population the anomaly is reported for, and a drift may be
contamination or composition. At 0.90 the set is fixed large fields, and any
drift is contamination alone.

What to read: if the anomaly is flat from the strict end down to the chosen
threshold, the choice is not a lever on the result — which is a stronger
statement than any claim that 0.9 is optimal. Where it starts to bend is where
contamination begins to matter.
"""

from __future__ import annotations

import logging

import matplotlib.pyplot as plt

from drought.analysis.sensitivity import load, sweep
from drought.core.config import CFG
from drought.viz.cuts import _save, _style

__all__ = ["figure"]

log = logging.getLogger(__name__)


def figure(site: str, regime: str = "intensive", min_px: int = 200):
    """The sensitivity figure for one site and system; returns the sweep table."""
    C = CFG
    cur = C.project.current_year
    thr_used = C.mask.purity_min
    t = sweep(load(site, regime), min_px=min_px)
    if t.empty:
        log.info("not enough pixels at any threshold")
        return None

    ref = t.iloc[(t.threshold - thr_used).abs().argmin()]
    strict = t.iloc[-1]

    plt.style.use("dark_background")
    fig, axes = plt.subplots(2, 1, figsize=(7.6, 7.4), dpi=180, sharex=True,
                             constrained_layout=True,
                             gridspec_kw={"height_ratios": [2.2, 1]})

    ax = axes[0]
    ax.plot(t.threshold, t.anomaly, "o-", ms=4, lw=2.2, color="#ff3b30")
    ax.axvline(thr_used, color="white", ls="--", lw=1.6,
               label=f"threshold in use ({thr_used:.2f})")
    ax.axhline(strict.anomaly, color="#8b949e", ls=":", lw=1.4,
               label=f"strictest available ({strict.threshold:.2f})")
    ax.fill_between(t.threshold, strict.anomaly * 0.95, strict.anomaly * 1.05,
                    color="#8b949e", alpha=0.18, lw=0,
                    label="within 5% of strictest")
    ax.set_ylabel(f"{cur} anomaly, {C.detect.index}", fontsize=10, color="white")
    ax.set_title("Does the answer depend on the threshold?", fontsize=11.5,
                 fontweight="bold", color="white", loc="left")
    ax.legend(fontsize=8, frameon=True, facecolor="#161b22", edgecolor="none",
              labelcolor="white", loc="lower right")
    ax.text(0.03, 0.97,
            f"at {ref.threshold:.2f}:  {ref.anomaly:+.4f}\n"
            f"at {strict.threshold:.2f}:  {strict.anomaly:+.4f}\n"
            f"difference:  {100 * abs(ref.anomaly - strict.anomaly) / abs(strict.anomaly):.1f}% "
            "of the anomaly",
            transform=ax.transAxes, ha="left", va="top", fontsize=9,
            color="#c9d1d9",
            bbox=dict(facecolor="#161b22", edgecolor="#30363d",
                      boxstyle="round,pad=0.4"))

    ax = axes[1]
    ax.plot(t.threshold, t.n_px, "o-", ms=4, lw=2.0, color="#58a6ff")
    ax.axvline(thr_used, color="white", ls="--", lw=1.6)
    ax.set_yscale("log")
    ax.set_xlabel("purity threshold", fontsize=10, color="white")
    ax.set_ylabel("cells contributing", fontsize=10, color="#58a6ff")
    ax.tick_params(axis="y", colors="#58a6ff")
    ax.set_title("How much data is left, and how close to the true area",
                 fontsize=11, fontweight="bold", color="white", loc="left")

    _style(axes)
    fig.suptitle(f"Purity threshold sensitivity  ·  {site}  ·  {regime}",
                 fontsize=11.5, fontweight="bold", color="#ff3b30")
    strict_p = C.sensitivity.strict_purity
    who = ("Every parcel enters (sensitivity.strict_purity 0): at high thresholds small fields "
           "drop out, so a drift can be contamination or composition."
           if not strict_p else
           f"Fixed parcel set (a cell >= {strict_p:g} required): the threshold changes "
           "which CELLS are used, never which FIELDS.")
    fig.supxlabel(f"{who}\nA flat curve means the choice is not a lever on the result.",
                  fontsize=7, color="#6e7681")
    _save(fig, site, f"purity_sensitivity_{regime}")
    plt.close(fig)

    lines = [f"{'threshold':>10}{'anomaly':>11}{'vs strictest':>14}{'cells':>10}"]
    for _, r in t.iterrows():
        dv = 100 * abs(r.anomaly - strict.anomaly) / abs(strict.anomaly)
        mark = "  <- in use" if abs(r.threshold - thr_used) < 1e-9 else ""
        lines.append(f"{r.threshold:>10.2f}{r.anomaly:>11.4f}{dv:>13.1f}%"
                     f"{int(r.n_px):>10,}{mark}")
    log.info("\n".join(lines))
    return t
