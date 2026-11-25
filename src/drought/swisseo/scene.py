"""One swissEO mosaic inside one site window: indices and quality layers.

The quality layers are kept as swissEO delivers them, and the rule that
decides which pixels are usable is applied later (`usable`), so it can change
without a new download:
  - CLOUDMASK-10M (OmniCloudMask): 0 clear, 1 thick cloud, 2 thin cloud,
    3 cloud shadow; usable if in swisseo.cloud_clear.
  - TERRAINMASK-10M: solar incidence angle 0-180 deg, 200 cast shadow, 255 no
    data; usable if not in swisseo.terrain_drop and below
    swisseo.terrain_max_incidence.
  - SCL-20M (scene classification): usable if not in swisseo.scl_drop
    (0 no data, 11 snow).
Under thick cloud and cloud shadow (NEVER_USABLE) no rule could use a pixel,
so the indices are not kept there.

The indices, from the same bands as before:
    NDVI  = (B8  - B4)  / (B8  + B4)          10 m
    kNDVI = tanh(NDVI^2)                       (sigma = 0.5 (NIR + Red))
    NDMI  = (B8A - B11) / (B8A + B11)          20 m, repeated 2 x 2
    NDRE  = (B8A - B5)  / (B8A + B5)           20 m, repeated 2 x 2
Reflectance = DN * swisseo.scale + swisseo.offset. The offset matters: a
normalised difference of DNs with the offset left in is a different number.
"""

from __future__ import annotations

import contextlib

import numpy as np

from ..config import CFG
from .grid import Window

__all__ = ["BANDS", "INDICES", "NEVER_USABLE", "QUALITY", "gdal_env", "indices",
           "quality", "read_band", "usable"]

INDICES = ("NDVI", "kNDVI", "NDMI", "NDRE")
BANDS = {"B4": "b04_10m", "B8": "b08_10m", "B5": "b05_20m", "B8A": "b8a_20m",
         "B11": "b11_20m"}
CLOUD, TERRAIN, SCL = "cloudmask_10m", "terrainmask_10m", "scl_20m"

# GDAL settings for reading COGs over HTTP: no directory listing, range
# requests merged, a block cache shared between the bands of a mosaic.
_GDAL = dict(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
             CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif",
             GDAL_HTTP_MERGE_CONSECUTIVE_RANGES="YES",
             GDAL_HTTP_MULTIPLEX="YES", GDAL_HTTP_VERSION="2",
             GDAL_HTTP_MAX_RETRY="4", GDAL_HTTP_RETRY_DELAY="2",
             VSI_CACHE="TRUE", GDAL_CACHEMAX=256)


@contextlib.contextmanager
def gdal_env():
    import rasterio
    with rasterio.Env(**_GDAL):
        yield


def read_band(href: str, win: Window, fill: int = 0) -> np.ndarray:
    """One band inside the window, on the window's 10 m grid.

    A 20 m band is read on its own grid (the 20 m cells covering the window)
    and each cell repeated 2 x 2, then cut to the window. Cells outside the
    mosaic come back as `fill`.
    """
    import rasterio
    from rasterio.windows import from_bounds

    with rasterio.open(href) as ds:
        res = ds.transform.a
        f = int(round(res / win.res))
        if f == 1:
            rw = from_bounds(*win.bounds, transform=ds.transform)
            rw = rw.round_offsets().round_lengths()
            return ds.read(1, window=rw, boundless=True, fill_value=fill)
        if f < 1 or abs(res - f * win.res) > 1e-6:
            raise ValueError(f"{href}: {res} m is not a multiple of {win.res} m")
        x0, y0, x1, y1 = win.bounds
        ox, oy = ds.transform.c, ds.transform.f
        gx0 = ox + np.floor((x0 - ox) / res) * res
        gy1 = oy - np.floor((oy - y1) / res) * res
        gx1 = ox + np.ceil((x1 - ox) / res) * res
        gy0 = oy - np.ceil((oy - y0) / res) * res
        rw = from_bounds(gx0, gy0, gx1, gy1, transform=ds.transform)
        rw = rw.round_offsets().round_lengths()
        a = ds.read(1, window=rw, boundless=True, fill_value=fill)
    a = np.repeat(np.repeat(a, f, axis=0), f, axis=1)
    r0 = int(round((gy1 - y1) / win.res))
    c0 = int(round((x0 - gx0) / win.res))
    return a[r0:r0 + win.height, c0:c0 + win.width]


NEVER_USABLE = (1, 3)        # thick cloud, cloud shadow
QUALITY = ("CLOUD", "TERRAIN", "SCL")


def quality(item, win: Window) -> dict[str, np.ndarray]:
    """The three quality layers as delivered, on the window grid (uint8)."""
    return {"CLOUD": read_band(item.href(CLOUD), win, fill=255),
            "TERRAIN": read_band(item.href(TERRAIN), win, fill=255),
            "SCL": read_band(item.href(SCL), win, fill=0)}


def usable(q: dict[str, np.ndarray]) -> np.ndarray:
    """Which pixels the current rules accept: clear, lit, not snow."""
    c = CFG.swisseo
    ok = np.isin(q["CLOUD"], list(c.cloud_clear))
    t = q["TERRAIN"]
    ok &= ~np.isin(t, list(c.terrain_drop)) & (t < float(c.terrain_max_incidence))
    ok &= ~np.isin(q["SCL"], list(c.scl_drop))
    return ok


def _refl(dn: np.ndarray) -> np.ndarray:
    c = CFG.swisseo
    r = dn.astype("float32") * float(c.scale) + float(c.offset)
    r[dn == 0] = np.nan                      # swissEO's no-data
    return r


def _nd(a, b):
    with np.errstate(invalid="ignore", divide="ignore"):
        return (a - b) / (a + b)


def indices(item, win: Window, keep: np.ndarray) -> dict[str, np.ndarray]:
    """NDVI, NDMI and NDRE on the window grid, NaN where not kept or not covered."""
    r = {k: _refl(read_band(item.href(v), win)) for k, v in BANDS.items()}
    out = {"NDVI": _nd(r["B8"], r["B4"]), "NDMI": _nd(r["B8A"], r["B11"]),
           "NDRE": _nd(r["B8A"], r["B5"])}
    for v in out.values():
        v[~keep] = np.nan
    return out
