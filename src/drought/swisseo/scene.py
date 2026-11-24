"""One swissEO mosaic inside one site window: the usable pixels and the indices.

Usable means all of:
  - covered: every band read has a reflectance (0 is swissEO's no-data);
  - clear: CLOUDMASK-10M (OmniCloudMask) in swisseo.cloud_clear;
  - lit: TERRAINMASK-10M not a cast shadow or no-data, and a solar incidence
    angle below swisseo.terrain_max_incidence;
  - not snow: SCL-20M not in swisseo.scl_drop.

The indices are those of the Earth Engine backend, from the same bands:
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

__all__ = ["BANDS", "INDICES", "clear_mask", "gdal_env", "indices", "read_band"]

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


def clear_mask(item, win: Window) -> np.ndarray:
    """Clear and lit and not snow, from the three quality layers. Cheap: uint8."""
    c = CFG.swisseo
    cloud = read_band(item.href(CLOUD), win, fill=255)
    ok = np.isin(cloud, list(c.cloud_clear))
    terr = read_band(item.href(TERRAIN), win, fill=255)
    ok &= ~np.isin(terr, list(c.terrain_drop)) & (terr < float(c.terrain_max_incidence))
    scl = read_band(item.href(SCL), win, fill=0)
    ok &= ~np.isin(scl, list(c.scl_drop))
    return ok


def _refl(dn: np.ndarray) -> np.ndarray:
    c = CFG.swisseo
    r = dn.astype("float32") * float(c.scale) + float(c.offset)
    r[dn == 0] = np.nan                      # swissEO's no-data
    return r


def _nd(a, b):
    with np.errstate(invalid="ignore", divide="ignore"):
        return (a - b) / (a + b)


def indices(item, win: Window, ok: np.ndarray) -> dict[str, np.ndarray]:
    """The four indices on the window grid, NaN wherever not usable."""
    r = {k: _refl(read_band(item.href(v), win)) for k, v in BANDS.items()}
    ndvi = _nd(r["B8"], r["B4"])
    out = {"NDVI": ndvi, "kNDVI": np.tanh(ndvi ** 2),
           "NDMI": _nd(r["B8A"], r["B11"]), "NDRE": _nd(r["B8A"], r["B5"])}
    for v in out.values():
        v[~ok] = np.nan
    return out
