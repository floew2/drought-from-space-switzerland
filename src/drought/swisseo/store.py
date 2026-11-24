"""The downloaded data, kept on disk: one GeoTIFF per mosaic and site.

Each file holds NDVI, NDMI and NDRE on the site window (10 m, LV95), as int16
scaled by 10 000, with -32768 as no data. A cell holds a value only if it is
usable (clear, lit, not snow: scene.clear_mask) and grassland: purity > 0 in
any management group of the mask. Every purity threshold therefore works from
these files without a new download; only a mask that adds grassland outside
the stored cells needs one (delete the site's folder).

    <paths.swisseo>/scenes/<site>/<mosaic id>.tif
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..config import series_root
from . import grid, scene

__all__ = ["BANDS", "fetch", "load", "path"]

BANDS = ("NDVI", "NDMI", "NDRE")
SCALE = 10_000
NODATA = -32768


def path(site: str, item_id: str) -> Path:
    return series_root() / "scenes" / site / f"{item_id}.tif"


def fetch(item, win: grid.Window, grass: np.ndarray, site: str) -> Path:
    """Download one mosaic for one site and store it. Returns the file."""
    import rasterio

    with scene.gdal_env():
        ok = scene.clear_mask(item, win) & grass
        ix = scene.indices(item, win, ok)
    arr = np.full((len(BANDS), *win.shape), NODATA, np.int16)
    for i, b in enumerate(BANDS):
        v = ix[b]
        good = ~np.isnan(v)
        arr[i][good] = np.round(np.clip(v[good], -1, 1) * SCALE).astype(np.int16)
    p = path(site, item.id)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".part.tif")
    with rasterio.open(tmp, "w", driver="GTiff", width=win.width, height=win.height,
                       count=len(BANDS), dtype="int16", nodata=NODATA, crs="EPSG:2056",
                       transform=win.transform, compress="deflate", predictor=2,
                       tiled=True, blockxsize=256, blockysize=256) as ds:
        ds.write(arr)
        ds.descriptions = BANDS
        ds.update_tags(datetime=item.datetime.isoformat(), source=item.href(scene.CLOUD))
    tmp.replace(p)
    return p


def load(p: Path, win: grid.Window) -> dict[str, np.ndarray]:
    """The stored indices as float32, NaN where no data, plus kNDVI."""
    import rasterio

    with rasterio.open(p) as ds:
        if (ds.height, ds.width) != win.shape or \
                tuple(round(v, 3) for v in ds.transform[:6]) != \
                tuple(round(v, 3) for v in win.transform[:6]):
            raise ValueError(f"{p}: stored on another window; delete the site's "
                             "scenes folder and run again")
        a = ds.read().astype("float32")
    a[a == NODATA] = np.nan
    out = {b: a[i] / SCALE for i, b in enumerate(BANDS)}
    out["kNDVI"] = np.tanh(out["NDVI"] ** 2)
    return out
