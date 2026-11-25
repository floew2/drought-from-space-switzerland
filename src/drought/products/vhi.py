"""swisstopo's Vegetation Health Index on our grassland cells: the cross-check.

VHI (swissEO VHI, with MeteoSwiss) combines a vegetation condition index from
NDVI and a temperature condition index from land-surface temperature, each
scaled between the 5th and 95th percentile of 1991-2020, weighted equally.
0-100, low is stressed. One map per day, each built from the last 7 days. Its
reference comes from Landsat (30 m) and its temperature from Meteosat, so its
10 m cells do not carry 10 m information throughout: it is read here as an
independent opinion, not as a measurement of the same thing as NDMI.

Read from the stored maps (`drought fetch`) on the same cells as the
condition series (purity >= mask.purity_min), current year, site window:

    <paths.series>/<site>/vhi_series.csv    per system and date: vhi_mean, n_px
    <paths.series>/<site>/vhi_parcels.csv   per parcel: mean of the daily parcel
                                            means, days observed
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from drought.core.config import CFG, series_dir, window
from drought.products.common import cells, parcels, rasterize, write
from drought.sources import store

__all__ = ["output", "run"]

log = logging.getLogger(__name__)


def output(site: str, what: str = "series") -> Path:
    return series_dir(site) / f"vhi_{what}.csv"


def run(site: str) -> list[Path]:
    """Compute the VHI series and per-parcel means of one site."""
    c = cells(site)
    p = parcels(site, c.win)
    year = int(CFG.project.current_year)
    labs, tabs = {}, {}
    for r in CFG.analysis.regimes:
        sub = p[p.lnf_code.isin(list(getattr(CFG.landuse, r)))].reset_index(drop=True)
        lab = rasterize(sub.geometry, np.arange(1, len(sub) + 1), c.win)
        lab[~c.regimes[r]] = 0
        labs[r], tabs[r] = lab, sub

    rows: list[list] = []
    per: dict[str, list[pd.Series]] = {r: [] for r in labs}
    for e, a in store.iterate(site, year, c.win, doy=window(site), kind="vhi"):
        for r, lab in labs.items():
            ok = (lab > 0) & ~np.isnan(a)
            if not ok.any():
                continue
            rows.append([site, r, f"{e.datetime:%Y-%m-%d}", e.doy, float(a[ok].mean()),
                         int(ok.sum())])
            per[r].append(pd.Series(a[ok]).groupby(lab[ok]).mean())
    series = pd.DataFrame(rows, columns=["site", "regime", "date", "doy", "vhi_mean", "n_px"])
    frames = []
    for r, ss in per.items():
        if not ss:
            continue
        d = pd.concat(ss, axis=1)
        t = tabs[r].iloc[d.index - 1]
        frames.append(pd.DataFrame({"parcel_id": t.parcel_id.values,
                                    "lnf_code": t.lnf_code.values, "regime": r,
                                    "vhi_mean": d.mean(axis=1).values,
                                    "n_days": d.notna().sum(axis=1).values}))
    parcels_df = (pd.concat(frames, ignore_index=True) if frames else
                  pd.DataFrame(columns=["parcel_id", "lnf_code", "regime", "vhi_mean", "n_days"]))
    return [write(series, output(site), "vhi"),
            write(parcels_df, output(site, "parcels"), "vhi")]
