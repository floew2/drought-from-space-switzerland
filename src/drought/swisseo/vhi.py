"""swisstopo's Vegetation Health Index on our grassland cells: the cross-check.

VHI (swissEO VHI, with MeteoSwiss) combines a vegetation condition index from
NDVI and a temperature condition index from land-surface temperature, each
scaled between the 5th and 95th percentile of 1991-2020, weighted equally.
0-100, low is stressed; 110 and 255 are no data. One map per day, each built
from the last 7 days. Its reference comes from Landsat (30 m) and its
temperature from Meteosat, so its 10 m cells do not carry 10 m information
throughout; it is read here as an independent opinion, not as a measurement
of the same thing as NDMI.

Read on the same cells as the trajectory (purity >= mask.purity_min), for the
current year, over each site's sensitivity window:

    vhi/vhi_<site>_<year>_series.csv    per management system and date: mean, n_px
    vhi/vhi_<site>_<year>_parcels.csv   per parcel: season mean of the daily
                                        parcel means, days observed
"""

from __future__ import annotations

import datetime as dt
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import CFG, ROOT, figures_dir, series_root
from . import extract, grid, scene, stac

__all__ = ["plot", "run"]


def _read(item: stac.Item, win: grid.Window) -> np.ndarray:
    a = scene.read_band(item.href(CFG.swisseo.vhi_asset), win, fill=255).astype("float32")
    a[np.isin(a, list(CFG.swisseo.vhi_nodata)) | (a > 100)] = np.nan
    return a


def run(sites=None, year: int | None = None, workers: int = 4, session=None,
        parcels=None) -> list[Path]:
    """Write the VHI series and per-parcel means; return the files."""
    from shapely.geometry import box

    year = int(year or CFG.project.current_year)
    sites = sites or list(CFG.gee.sensitivity.windows.to_dict())
    p_all = extract.parcels_table() if parcels is None else parcels
    out_dir = series_root() / "vhi"
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for site in sites:
        win = grid.site_window(site)
        m = grid.read_mask(win)
        pmin = float(CFG.mask.purity_min)
        w = getattr(CFG.gee.sensitivity.windows, site)
        start = dt.date(year, 1, 1) + dt.timedelta(days=int(w[0]) - 1)
        end = dt.date(year, 1, 1) + dt.timedelta(days=int(w[1]) - 1)
        items = stac.search(getattr(CFG.sites, site).bbox_wgs84, start, end,
                            collection=CFG.swisseo.vhi_collection, session=session)
        p = p_all[p_all.intersects(box(*win.bounds))].reset_index(drop=True)
        labs, tabs = {}, {}
        for r in CFG.analysis.regimes:
            sub = p[p.lnf_code.isin(list(getattr(CFG.landuse, r)))].reset_index(drop=True)
            lab = extract._rasterize(sub.geometry, np.arange(1, len(sub) + 1), win)
            lab[~(m[r] >= pmin)] = 0
            labs[r], tabs[r] = lab, sub

        def one(it, win=win):
            with scene.gdal_env():
                return it, _read(it, win)

        with ThreadPoolExecutor(max_workers=max(1, workers)) as ex:
            res = list(ex.map(one, items))

        rows, per = [], {r: [] for r in labs}
        for it, a in res:
            for r, lab in labs.items():
                cells = (lab > 0) & ~np.isnan(a)
                if not cells.any():
                    continue
                rows.append([site, r, f"{it.datetime:%Y-%m-%d}", it.doy,
                             float(a[cells].mean()), int(cells.sum())])
                s = pd.Series(a[cells]).groupby(lab[cells]).mean()
                per[r].append(s)
        series = pd.DataFrame(rows, columns=["site", "regime", "date", "doy", "vhi_mean", "n_px"])
        f1 = out_dir / f"vhi_{site}_{year}_series.csv"
        series.to_csv(f1, index=False)
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
        f2 = out_dir / f"vhi_{site}_{year}_parcels.csv"
        (pd.concat(frames, ignore_index=True) if frames else
         pd.DataFrame(columns=["parcel_id", "lnf_code", "regime", "vhi_mean", "n_days"])
         ).to_csv(f2, index=False)
        for f in (f1, f2):
            print(f"  -> {f.relative_to(ROOT) if f.is_relative_to(ROOT) else f}")
        written += [f1, f2]
    return written


def plot(year: int | None = None, sites=None, save: bool = True):
    """VHI on our pure grassland cells through the season, one panel per site."""
    import matplotlib.pyplot as plt

    year = int(year or CFG.project.current_year)
    sites = sites or list(CFG.gee.sensitivity.windows.to_dict())
    fig, axes = plt.subplots(1, len(sites), figsize=(4.2 * len(sites), 3.4),
                             sharey=True, squeeze=False)
    for ax, site in zip(axes[0], sites):
        f = series_root() / "vhi" / f"vhi_{site}_{year}_series.csv"
        ax.axhspan(40, 50, color="#cbffca", alpha=0.6, lw=0, label="normal (40-50)")
        ax.axhspan(0, 30, color="#f5cd85", alpha=0.35, lw=0, label="stressed (< 30)")
        if f.exists():
            d = pd.read_csv(f, parse_dates=["date"])
            for r, g in d.groupby("regime"):
                ax.plot(g.date, g.vhi_mean, lw=1.4, label=r)
        ax.set_title(getattr(CFG.sites, site).label, fontsize=10)
        ax.set_ylim(0, 100)
        ax.tick_params(axis="x", labelrotation=30, labelsize=8)
    axes[0][0].set_ylabel("VHI (swissEO)")
    axes[0][-1].legend(fontsize=7, loc="lower right")
    fig.suptitle(f"Vegetation Health Index on pure grassland cells, {year}", fontsize=11)
    fig.tight_layout()
    if save:
        out = figures_dir() / f"vhi_{year}.png"
        fig.savefig(out, dpi=150)
        print(f"  -> {out.relative_to(ROOT) if out.is_relative_to(ROOT) else out}")
    return fig
