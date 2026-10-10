"""Reference data: saving, loading and validating hand labels.

Labels are expensive to produce and must outlive the kernel. Once saved, every
later change to the detector is scored against the SAME ground truth, so an
improvement is measured rather than asserted.

One file per site. A detector validated on Linthebene has not been validated on
Rheintal, and pooling the labels would hide that.
"""

from __future__ import annotations

import json

import pandas as pd

from .config import CFG, labels_path


def save_labels(labels: dict, site: str) -> None:
    p = labels_path(site)
    with open(p, "w") as f:
        json.dump({k: sorted(int(x) for x in v) for k, v in labels.items()},
                  f, indent=1, sort_keys=True)
    n = sum(len(v) for v in labels.values())
    print(f"wrote {len(labels)} parcel-years, {n} cuts -> {p}")


def load_labels(site: str) -> dict:
    p = labels_path(site)
    if not p.exists():
        print(f"no labels yet for {site} ({p})")
        return {}
    with open(p) as f:
        labels = json.load(f)
    n = sum(len(v) for v in labels.values())
    print(f"loaded {len(labels)} parcel-years, {n} cuts from {p}")
    return labels


def check_labels(labels: dict, df: pd.DataFrame, cfg=None) -> dict:
    """Flag labels the detector could never match, before they distort scoring.

    A label outside the season, or on a parcel-year absent from the data, counts
    as a permanent miss and silently depresses recall no matter what the
    detector does. Better to see them than to tune against them.
    """
    C = cfg or CFG
    lo, hi = C.detect.season
    ok, dropped = {}, []
    for key, v in labels.items():
        pid, yr = key.split("|")
        g = df[(df.parcel_id == pid) & (df.year == int(yr))]
        if g.empty:
            dropped.append((key, "parcel-year not in data")); continue
        keep = [d for d in v if lo <= d <= hi]
        if len(keep) < len(v):
            dropped.append((key, f"outside season ({lo}, {hi}): "
                                 f"{[d for d in v if d not in keep]}"))
        if keep:
            ok[key] = keep
    for key, why in dropped:
        print(f"  [!] {key}: {why}")
    print(f"  {len(ok)} parcel-years usable, "
          f"{sum(len(v) for v in ok.values())} cuts")
    yrs = pd.Series([int(k.split("|")[1]) for k in ok]).value_counts().sort_index()
    print("  parcel-years per year: "
          + ", ".join(f"{y}:{n}" for y, n in yrs.items()))
    return ok
