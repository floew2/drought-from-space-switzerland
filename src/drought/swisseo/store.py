"""The downloaded data, kept on disk: one GeoTIFF per mosaic and site.

Each file holds, on the site window (10 m, LV95), int16:
    NDVI, NDMI, NDRE      scaled by 10 000; -32768 where not covered, or under
                          thick cloud or cloud shadow (scene.NEVER_USABLE)
    CLOUD, TERRAIN, SCL   swissEO's quality layers, unchanged
Nothing in a file depends on the mask or on the cloud rule: `load` applies the
current rule (scene.usable), so any mask, purity threshold, parcel set or
cloud rule works from the same files without a new download.

    <paths.scenes>/<site>/<mosaic id>.tif
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np

from ..config import series_root
from . import grid, scene

__all__ = ["BANDS", "fetch", "load", "path"]

INDEX_BANDS = ("NDVI", "NDMI", "NDRE")
BANDS = (*INDEX_BANDS, *scene.QUALITY)
SCALE = 10_000
NODATA = -32768
ATTEMPTS = 4          # a read cut off in transit is retried; waits 5, 20, 45 s


def path(site: str, item_id: str) -> Path:
    return series_root() / "scenes" / site / f"{item_id}.tif"


def _download(item, win: grid.Window) -> tuple[dict, dict]:
    from rasterio.errors import RasterioError

    for attempt in range(1, ATTEMPTS + 1):
        try:
            with scene.gdal_env():
                q = scene.quality(item, win)
                keep = ~np.isin(q["CLOUD"], scene.NEVER_USABLE)
                return scene.indices(item, win, keep), q
        except (RasterioError, OSError):
            if attempt == ATTEMPTS:
                raise
            time.sleep(5 * attempt ** 2)
    raise AssertionError("unreachable")


def fetch(item, win: grid.Window, site: str) -> Path:
    """Download one mosaic for one site and store it. Returns the file.

    Raises:
        rasterio.errors.RasterioError: The mosaic could not be read after
            ATTEMPTS tries.
    """
    import rasterio

    ix, q = _download(item, win)
    arr = np.full((len(BANDS), *win.shape), NODATA, np.int16)
    for i, b in enumerate(INDEX_BANDS):
        v = ix[b]
        good = ~np.isnan(v)
        arr[i][good] = np.round(np.clip(v[good], -1, 1) * SCALE).astype(np.int16)
    for i, b in enumerate(scene.QUALITY, start=len(INDEX_BANDS)):
        arr[i] = q[b].astype(np.int16)
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
    """The stored indices as float32, NaN wherever the current rule rejects the
    pixel (scene.usable) or there is no data; plus kNDVI."""
    import rasterio

    with rasterio.open(p) as ds:
        if (ds.height, ds.width) != win.shape or \
                tuple(round(v, 3) for v in ds.transform[:6]) != \
                tuple(round(v, 3) for v in win.transform[:6]):
            raise ValueError(f"{p}: stored on another window; delete the site's "
                             "scenes folder and run again")
        if tuple(ds.descriptions) != BANDS:
            raise ValueError(f"{p}: bands {ds.descriptions}, expected {BANDS}: stored by "
                             "an older version; delete the scenes folder and run again")
        a = ds.read()
    q = {b: a[i].astype(np.uint8) for i, b in enumerate(scene.QUALITY, start=len(INDEX_BANDS))}
    ok = scene.usable(q)
    out = {}
    for i, b in enumerate(INDEX_BANDS):
        v = a[i].astype("float32")
        v[(a[i] == NODATA) | ~ok] = np.nan
        out[b] = v / SCALE
    out["kNDVI"] = np.tanh(out["NDVI"] ** 2)
    return out
