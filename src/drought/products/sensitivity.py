"""How the anomaly depends on the purity threshold.

Per management system and year: the per-cell median of sensitivity.index over
the site's window, then its mean per purity class (n_classes bins of 0 to 1).
The cells are those inside parcels of the system with purity > 0, in parcels
holding a cell at sensitivity.strict_purity (0.0 = every parcel). A mosaic
whose share of unusable grassland exceeds sensitivity.max_site_cloud is left
out of the composite.

    <paths.series>/<site>/sensitivity_<regime>.csv
        year, pclass, purity_lo, mean, n_px
"""

from __future__ import annotations

import logging
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from drought.core.config import CFG, series_dir, window, years
from drought.products.common import cells, parcels, rasterize, write
from drought.sources import store

__all__ = ["COLUMNS", "output", "run"]

log = logging.getLogger(__name__)
COLUMNS = ["year", "pclass", "purity_lo", "mean", "n_px"]


def output(site: str, regime: str) -> Path:
    return series_dir(site) / f"sensitivity_{regime}.csv"


def run(site: str) -> list[Path]:
    """Compute the sensitivity series of one site, one file per system."""
    c = cells(site)
    p = parcels(site, c.win)
    s = CFG.sensitivity
    n_cls, strict = int(s.n_classes), float(s.strict_purity)
    idx, cls = {}, {}
    for r in CFG.analysis.regimes:
        sub = p[p.lnf_code.isin(list(getattr(CFG.landuse, r)))]
        lab = rasterize(sub.geometry, np.arange(1, len(sub) + 1), c.win)
        if strict > 0:
            hit = np.unique(lab[(lab > 0) & (c.purity[r] >= strict)])
            lab = np.where(np.isin(lab, hit), lab, 0)
        flat = np.flatnonzero((lab > 0) & (c.purity[r] > 0))
        idx[r] = flat
        cls[r] = np.minimum(np.floor(c.purity[r].ravel()[flat] * n_cls), n_cls - 1).astype(int)

    frames: dict[str, list[pd.DataFrame]] = {r: [] for r in idx}
    for year in years():
        stack: dict[str, list[np.ndarray]] = {r: [] for r in idx}
        for _, ix in store.iterate(site, year, c.win, doy=window(site)):
            if c.cloud_pct(ix) > float(s.max_site_cloud):
                continue
            v = ix[s.index].ravel()
            for r, flat in idx.items():
                stack[r].append(v[flat].astype("float32"))
        for r, st in stack.items():
            if not st:
                continue
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)       # all-NaN cells
                comp = np.nanmedian(np.vstack(st), axis=0)
            d = pd.DataFrame({"pclass": cls[r], "v": comp}).dropna()
            g = d.groupby("pclass")["v"].agg(["mean", "count"]).reset_index()
            frames[r].append(pd.DataFrame({"year": year, "pclass": g.pclass,
                                           "purity_lo": g.pclass / n_cls,
                                           "mean": g["mean"], "n_px": g["count"]}))
    out = []
    for r, fr in frames.items():
        if fr:
            out.append(write(pd.concat(fr, ignore_index=True)[COLUMNS], output(site, r),
                             "sensitivity"))
        else:
            log.warning("  [!] %s/%s: no usable mosaic in the window", site, r)
    return out
