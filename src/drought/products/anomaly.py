"""Per-field anomaly: each parcel's current year against its own baseline.

For every parcel of a management system inside the AOI, on its pure cells
(purity >= mask.purity_min for that system):
  1. per stored mosaic in the site's window, the parcel mean of anomaly.index,
     if at least half of its pure cells are usable (a partly clouded parcel
     samples a different part of it);
  2. per year, the mean of those, if there are at least anomaly.min_obs;
  3. z = (current year - median of the baseline years) / their standard
     deviation, if at least anomaly.min_baseline_years baseline years remain.

    <paths.series>/<site>/anomaly.csv
        parcel_id, lnf_code, regime, n_cells, <year>..., n_baseline, baseline_median,
        baseline_sd, z
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from drought.core.config import CFG, series_dir, window, years
from drought.products.common import cells, parcels, rasterize, write
from drought.sources import store

__all__ = ["MIN_SHARE", "output", "run"]

log = logging.getLogger(__name__)
MIN_SHARE = 0.5          # share of a parcel's pure cells that must be usable in a mosaic


def output(site: str) -> Path:
    return series_dir(site) / "anomaly.csv"


def run(site: str) -> Path:
    """Compute the per-field anomaly of one site."""
    c = cells(site)
    p = parcels(site, c.win)
    a = CFG.anomaly
    cur = int(CFG.project.current_year)
    base_years = [int(y) for y in CFG.project.baseline_years]

    labels, tables = [], []
    offset = 0
    lab_all = np.zeros(c.win.shape, np.int32)
    for r in CFG.analysis.regimes:
        sub = p[p.lnf_code.isin(list(getattr(CFG.landuse, r)))].reset_index(drop=True)
        lab = rasterize(sub.geometry, np.arange(1, len(sub) + 1) + offset, c.win)
        lab[~c.regimes[r]] = 0
        lab_all = np.where(lab > 0, lab, lab_all)
        tables.append(sub.assign(regime=r, key=np.arange(1, len(sub) + 1) + offset))
        offset += len(sub)
        labels.append(lab)
    t = pd.concat(tables, ignore_index=True)
    flat = lab_all.ravel()
    has = flat > 0
    keys = flat[has]
    n_cells = np.bincount(keys, minlength=offset + 1)
    t["n_cells"] = n_cells[t.key]
    t = t[t.n_cells > 0].reset_index(drop=True)

    per_year = {}
    for year in years():
        sums = np.zeros(offset + 1)
        nobs = np.zeros(offset + 1)
        for _, ix in store.iterate(site, year, c.win, doy=window(site)):
            v = ix[a.index].ravel()[has]
            ok = ~np.isnan(v)
            cnt = np.bincount(keys[ok], minlength=offset + 1)
            tot = np.bincount(keys[ok], weights=v[ok], minlength=offset + 1)
            good = cnt >= MIN_SHARE * np.maximum(n_cells, 1)
            good &= n_cells > 0
            sums[good] += tot[good] / cnt[good]
            nobs[good] += 1
        with np.errstate(invalid="ignore", divide="ignore"):
            val = np.where(nobs >= int(a.min_obs), sums / nobs, np.nan)
        per_year[year] = val[t.key]

    vals = pd.DataFrame(per_year)
    base = vals[base_years]
    t["n_baseline"] = base.notna().sum(axis=1).to_numpy()
    t["baseline_median"] = base.median(axis=1).to_numpy()
    t["baseline_sd"] = base.std(axis=1, ddof=1).to_numpy()
    enough = (t.n_baseline >= int(a.min_baseline_years)) & vals[cur].notna().to_numpy()
    with np.errstate(invalid="ignore", divide="ignore"):
        z = (vals[cur].to_numpy() - t.baseline_median) / t.baseline_sd
    t["z"] = np.where(enough & (t.baseline_sd > 0), z, np.nan)
    out = pd.concat([t[["parcel_id", "lnf_code", "regime", "n_cells"]],
                     vals.rename(columns=str),
                     t[["n_baseline", "baseline_median", "baseline_sd", "z"]]], axis=1)
    log.info("  %s: %d parcels, %d with a z-score, median z %.2f", site, len(out),
             int(out.z.notna().sum()), float(np.nanmedian(out.z)) if out.z.notna().any()
             else float("nan"))
    return write(out, output(site), "anomaly")
