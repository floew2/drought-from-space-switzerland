"""The downloaded data, kept on disk, and the catalogue of what was downloaded.

Sentinel-2 (S2), one GeoTIFF per mosaic and site, on the site window
(10 m, LV95), int16:
    NDVI, NDMI, NDRE      scaled by 10 000; -32768 where not covered, or under
                          thick cloud or cloud shadow (scene.NEVER_USABLE)
    CLOUD, TERRAIN, SCL   swissEO's quality layers, unchanged
Nothing in a file depends on the mask or on the cloud rule: `load` applies the
current rule (scene.usable), so any mask, purity threshold, parcel set or cloud
rule works from the same files without a new download.

VHI, one single-band uint8 GeoTIFF per day and site, as delivered (0-100;
swisseo.vhi.nodata where there is none).

    <paths.scenes>/<site>/s2/<mosaic id>.tif
    <paths.scenes>/<site>/vhi/<map id>.tif
    <paths.scenes>/<site>/catalogue_<kind>_<year>.json   what the catalogue held
                                                          when it was fetched

The catalogue file is what every product iterates over, so computing needs no
network: it lists every mosaic of the year, and `complete` says whether all of
them are on disk.
"""

from __future__ import annotations

import datetime as dt
import json
import time
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from drought.core import grid
from drought.core.config import scenes_dir
from drought.sources import scene

__all__ = ["BANDS", "Entry", "catalogue", "complete", "fetch", "fetch_vhi", "iterate",
           "load", "load_vhi", "path", "write_catalogue"]

INDEX_BANDS = ("NDVI", "NDMI", "NDRE")
BANDS = (*INDEX_BANDS, *scene.QUALITY)
SCALE = 10_000
NODATA = -32768
ATTEMPTS = 4          # a read cut off in transit is retried; waits 5, 20, 45 s
KINDS = ("s2", "vhi")


@dataclass(frozen=True)
class Entry:
    """One mosaic in the catalogue file."""
    id: str
    datetime: dt.datetime

    @property
    def doy(self) -> int:
        return self.datetime.timetuple().tm_yday


def path(site: str, item_id: str, kind: str = "s2") -> Path:
    return scenes_dir(site) / kind / f"{item_id}.tif"


def _catalogue_file(site: str, year: int, kind: str) -> Path:
    return scenes_dir(site) / f"catalogue_{kind}_{year}.json"


def write_catalogue(site: str, year: int, items, kind: str = "s2") -> Path:
    """Record which mosaics the catalogue held for a site and year."""
    f = _catalogue_file(site, year, kind)
    f.write_text(json.dumps([{"id": i.id, "datetime": i.datetime.isoformat()} for i in items],
                            indent=0))
    return f


def catalogue(site: str, year: int, kind: str = "s2") -> list[Entry]:
    """The mosaics of a site and year, as fetched.

    Raises:
        FileNotFoundError: If that year was never fetched for the site.
    """
    f = _catalogue_file(site, year, kind)
    if not f.exists():
        raise FileNotFoundError(f"{site} {year}: no {kind} catalogue in {f.parent}; "
                                f"run `drought fetch --aoi {site}`")
    return [Entry(e["id"], dt.datetime.fromisoformat(e["datetime"]))
            for e in json.loads(f.read_text())]


def complete(site: str, year: int, kind: str = "s2") -> bool:
    """Every mosaic listed for that year is on disk."""
    try:
        return all(path(site, e.id, kind).exists() for e in catalogue(site, year, kind))
    except FileNotFoundError:
        return False


# --------------------------------------------------------------------------- #
# Download
# --------------------------------------------------------------------------- #
def _retry(read):
    from rasterio.errors import RasterioError

    for attempt in range(1, ATTEMPTS + 1):
        try:
            with scene.gdal_env():
                return read()
        except (RasterioError, OSError):
            if attempt == ATTEMPTS:
                raise
            time.sleep(5 * attempt ** 2)
    raise AssertionError("unreachable")


def _write(p: Path, arr: np.ndarray, win: grid.Window, names, nodata, tags: dict) -> Path:
    import rasterio

    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".part.tif")
    with rasterio.open(tmp, "w", driver="GTiff", width=win.width, height=win.height,
                       count=arr.shape[0], dtype=arr.dtype.name, nodata=nodata,
                       crs="EPSG:2056", transform=win.transform, compress="deflate",
                       predictor=2, tiled=True, blockxsize=256, blockysize=256) as ds:
        ds.write(arr)
        ds.descriptions = tuple(names)
        ds.update_tags(**tags)
    tmp.replace(p)
    return p


def fetch(item, win: grid.Window, site: str) -> Path:
    """Download one S2 mosaic for one site and store it. Returns the file.

    Raises:
        rasterio.errors.RasterioError: The mosaic could not be read after
            ATTEMPTS tries.
    """
    def read():
        q = scene.quality(item, win)
        keep = ~np.isin(q["CLOUD"], scene.NEVER_USABLE)
        return scene.indices(item, win, keep), q

    ix, q = _retry(read)
    arr = np.full((len(BANDS), *win.shape), NODATA, np.int16)
    for i, b in enumerate(INDEX_BANDS):
        v = ix[b]
        good = ~np.isnan(v)
        arr[i][good] = np.round(np.clip(v[good], -1, 1) * SCALE).astype(np.int16)
    for i, b in enumerate(scene.QUALITY, start=len(INDEX_BANDS)):
        arr[i] = q[b].astype(np.int16)
    return _write(path(site, item.id), arr, win, BANDS, NODATA,
                  {"datetime": item.datetime.isoformat(), "source": item.href(scene.CLOUD)})


def fetch_vhi(item, win: grid.Window, site: str) -> Path:
    """Download one VHI map for one site and store it as delivered."""
    from drought.core.config import CFG

    href = item.href(CFG.swisseo.vhi.asset)
    a = _retry(lambda: scene.read_band(href, win, fill=255)).astype(np.uint8)
    return _write(path(site, item.id, "vhi"), a[None], win, ("VHI",), 255,
                  {"datetime": item.datetime.isoformat(), "source": href})


# --------------------------------------------------------------------------- #
# Reading back
# --------------------------------------------------------------------------- #
def _read(p: Path, win: grid.Window, bands) -> np.ndarray:
    import rasterio

    with rasterio.open(p) as ds:
        if (ds.height, ds.width) != win.shape or \
                tuple(round(v, 3) for v in ds.transform[:6]) != \
                tuple(round(v, 3) for v in win.transform[:6]):
            raise ValueError(f"{p}: stored on another window; delete the site's "
                             "scenes folder and fetch again")
        if tuple(ds.descriptions) != tuple(bands):
            raise ValueError(f"{p}: bands {ds.descriptions}, expected {tuple(bands)}: "
                             "stored by an older version; delete the scenes folder and "
                             "fetch again")
        return ds.read()


def load(p: Path, win: grid.Window) -> dict[str, np.ndarray]:
    """The stored indices as float32, NaN wherever the current rule rejects the
    pixel (scene.usable) or there is no data; plus kNDVI."""
    a = _read(p, win, BANDS)
    q = {b: a[i].astype(np.uint8)
         for i, b in enumerate(scene.QUALITY, start=len(INDEX_BANDS))}
    ok = scene.usable(q)
    out = {}
    for i, b in enumerate(INDEX_BANDS):
        v = a[i].astype("float32")
        v[(a[i] == NODATA) | ~ok] = np.nan
        out[b] = v / SCALE
    out["kNDVI"] = np.tanh(out["NDVI"] ** 2)
    return out


def load_vhi(p: Path, win: grid.Window) -> np.ndarray:
    """VHI 0-100 as float32, NaN where swissEO has none."""
    from drought.core.config import CFG

    a = _read(p, win, ("VHI",))[0].astype("float32")
    a[np.isin(a, list(CFG.swisseo.vhi.nodata)) | (a > 100)] = np.nan
    return a


def iterate(site: str, year: int, win: grid.Window, doy: tuple[int, int] | None = None,
            kind: str = "s2") -> Iterator[tuple[Entry, Any]]:
    """The stored mosaics of a site and year, in time order, read from disk.

    Raises:
        FileNotFoundError: If the year was not fetched, or a mosaic it lists is
            missing (the fetch did not finish).
    """
    for e in catalogue(site, year, kind):
        if doy and not (doy[0] <= e.doy <= doy[1]):
            continue
        p = path(site, e.id, kind)
        if not p.exists():
            raise FileNotFoundError(f"{site} {year}: {e.id} not downloaded; "
                                    f"run `drought fetch --aoi {site}` to finish")
        yield e, (load(p, win) if kind == "s2" else load_vhi(p, win))
