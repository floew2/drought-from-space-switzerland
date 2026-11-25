"""Grassland condition: site-mean index series per management system.

For every stored mosaic in condition.doy: the mean of NDVI, kNDVI, NDMI and
NDRE over the cells whose purity reaches mask.purity_min, and how many of them
were usable. A date is kept if that count reaches condition.min_pixels.

    <paths.series>/<site>/condition.csv
        site, regime, year, date, doy, NDVI, kNDVI, NDMI, NDRE, n_px
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from drought.core.config import CFG, series_dir, years
from drought.products.common import cells, write
from drought.sources import scene, store

__all__ = ["COLUMNS", "output", "run"]

log = logging.getLogger(__name__)
COLUMNS = ["site", "regime", "year", "date", "doy", *scene.INDICES, "n_px"]


def output(site: str) -> Path:
    return series_dir(site) / "condition.csv"


def run(site: str) -> Path:
    """Compute the condition series of one site from its stored mosaics."""
    c = cells(site)
    lo, hi = (int(v) for v in CFG.condition.doy)
    floor = int(CFG.condition.min_pixels)
    rows = []
    for year in years():
        n_dates = 0
        for e, ix in store.iterate(site, year, c.win, doy=(lo, hi)):
            valid = ~np.isnan(ix["NDVI"])
            for regime, pure in c.regimes.items():
                n = int((pure & valid).sum())
                if n < floor:
                    continue
                with np.errstate(all="ignore"):
                    means = [float(np.nanmean(ix[i][pure])) for i in scene.INDICES]
                if not any(np.isnan(means)):
                    rows.append([site, regime, year, f"{e.datetime:%Y-%m-%d}", e.doy, *means, n])
                    n_dates += 1
        log.info("  %s %d: %d site-dates with enough clear grassland", site, year, n_dates)
    return write(pd.DataFrame(rows, columns=COLUMNS), output(site), "condition")
