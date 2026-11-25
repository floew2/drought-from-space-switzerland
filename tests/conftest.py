"""Shared fixtures: a tiny swissEO world on disk, served by a fake catalogue.

The world is one site, 40 x 40 cells of 10 m: a mask, two parcels, one mosaic
per year on 1 July (DOY 182, inside every window) and a half-cloudy one on
11 July of the current year. NDMI falls in the current year, so every product
has a drought to find.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pytest
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import box

X0, Y0, N = 2_720_000.0, 1_230_000.0, 40
BASE_DN = {"b04_10m": 1100, "b08_10m": 3100, "b05_20m": 1600, "b8a_20m": 3100, "b11_20m": 2100}


def refl(dn):
    return dn * 1e-4 - 0.1


def ndmi(b8a, b11):
    return (refl(b8a) - refl(b11)) / (refl(b8a) + refl(b11))


def tif(path, arr, x0, y0, res, dtype, desc=None):
    arr = np.asarray(arr)
    arr = arr[None] if arr.ndim == 2 else arr
    with rasterio.open(path, "w", driver="GTiff", width=arr.shape[2], height=arr.shape[1],
                       count=arr.shape[0], dtype=dtype, crs="EPSG:2056",
                       transform=from_origin(x0, y0, res, res)) as ds:
        ds.write(arr.astype(dtype))
        for i, d in enumerate(desc or [], 1):
            ds.set_band_description(i, d)
    return str(path)


def mosaic(folder, b11=2100, cloudy_rows=0, cloud_class=1):
    """One mosaic around the site: 10 m bands on a 10 m grid, 20 m bands on a
    20 m grid whose origin lies OFF the window's 20 m phase."""
    folder.mkdir(parents=True, exist_ok=True)
    gx, gy, n10 = X0 - 200, Y0 + 200, N + 40
    dn = {**BASE_DN, "b11_20m": b11}
    hrefs = {}
    for k, v in dn.items():
        if k.endswith("10m"):
            hrefs[k] = tif(folder / f"{k}.tif", np.full((n10, n10), v), gx, gy, 10, "uint16")
        else:
            hrefs[k] = tif(folder / f"{k}.tif", np.full((n10 // 2 + 2, n10 // 2 + 2), v),
                           gx - 10, gy + 10, 20, "uint16")
    hrefs["b03_10m"] = tif(folder / "b03.tif", np.full((n10, n10), 1500), gx, gy, 10, "uint16")
    cloud = np.zeros((n10, n10))
    cloud[:cloudy_rows] = cloud_class
    hrefs["cloudmask_10m"] = tif(folder / "cloud.tif", cloud, gx, gy, 10, "uint8")
    hrefs["terrainmask_10m"] = tif(folder / "terrain.tif", np.full((n10, n10), 40),
                                   gx, gy, 10, "uint8")
    hrefs["scl_20m"] = tif(folder / "scl.tif", np.full((n10 // 2 + 2, n10 // 2 + 2), 4),
                           gx - 10, gy + 10, 20, "uint8")
    return hrefs


class _Resp:
    def __init__(self, page): self.page = page
    def raise_for_status(self): pass
    def json(self): return self.page


class FakeCatalogue:
    """Answers items queries per collection with the mosaics in the date range."""

    def __init__(self):
        self.items = {}                     # collection -> {id: (iso date, hrefs)}
        self.calls = 0

    def add(self, collection, item_id, date, hrefs):
        self.items.setdefault(collection, {})[item_id] = (date, hrefs)

    def get(self, url, params=None, timeout=None):
        self.calls += 1
        coll = url.rstrip("/").split("/")[-2]
        a, b = (x[:10] for x in params["datetime"].split("/"))
        feats = [{"id": i, "properties": {"datetime": f"{d}T10:30:00Z"},
                  "assets": {f"x_mosaic_{d}t103000_{k}.tif": {"href": h}
                             for k, h in hrefs.items()}}
                 for i, (d, hrefs) in self.items.get(coll, {}).items() if a <= d <= b]
        return _Resp({"features": feats, "links": []})


@pytest.fixture
def world(tmp_path, monkeypatch):
    """Paths in tmp_path, a mask, a cadastre with ids, labels, and a catalogue."""
    from pyproj import Transformer

    from drought.core import config
    from drought.core.config import CFG
    from drought.core.provenance import record, write_sidecar
    from drought.sources import stac

    for k in ("scenes", "series", "masks", "labels", "figures", "tables", "parcel_ids"):
        monkeypatch.setattr(CFG.paths, k, str(tmp_path / k))
    monkeypatch.setattr(CFG.condition, "min_pixels", 10)
    monkeypatch.setattr(CFG.cuts, "min_area_ha", 0.5)
    monkeypatch.setattr(CFG.anomaly, "min_obs", 1)

    purity = np.full((4, N, N), 0.95, "float32")
    purity[:, :, :5] = 0.5                                     # an impure strip
    tif(config.mask_path(), purity, X0, Y0, 10, "float32",
        ["mown", "grazed", "intensive", "extensive"])
    write_sidecar(config.mask_path(), record("mask", config.mask_path()))
    bb = Transformer.from_crs(2056, 4326, always_xy=True).transform_bounds(
        X0 + 25, Y0 - N * 10 + 25, X0 + N * 10 - 25, Y0 - 25)
    monkeypatch.setattr(CFG.sites.linthebene, "bbox_wgs84", list(bb))

    parcels = gpd.GeoDataFrame(
        {"parcel_id": ["2720120_1229880_14400", "2720280_1229720_14400"],
         "lnf_code": [613, 611]},
        geometry=[box(X0 + 60, Y0 - 180, X0 + 180, Y0 - 60),
                  box(X0 + 220, Y0 - 340, X0 + 340, Y0 - 220)], crs=2056)
    from drought.sources import parcels as cadastre
    cadastre.path().parent.mkdir(parents=True, exist_ok=True)
    parcels.to_file(cadastre.path())

    cat = FakeCatalogue()
    cur = int(CFG.project.current_year)
    for y in config.years():
        b11 = 2600 if y == cur else 2100 + 20 * (y - 2019)      # drier in the current year
        cat.add(CFG.swisseo.collection, f"{y}-07-01t103000", f"{y}-07-01",
                mosaic(tmp_path / "img" / str(y), b11=b11))
    cat.add(CFG.swisseo.collection, f"{cur}-07-11t103000", f"{cur}-07-11",
            mosaic(tmp_path / "img" / "cloudy", b11=2600, cloudy_rows=40))
    vhi = np.full((N + 40, N + 40), 25, "uint8")
    vhi[:40] = 110                                              # no data, top 20 mask rows
    cat.add(CFG.swisseo.vhi.collection, f"{cur}-07-01t235959", f"{cur}-07-01",
            {CFG.swisseo.vhi.asset: tif(tmp_path / "vhi.tif", vhi, X0 - 200, Y0 + 200, 10,
                                        "uint8")})
    real_search = stac.search
    monkeypatch.setattr(stac, "search",
                        lambda *a, session=None, **k: real_search(*a, session=cat, **k))
    monkeypatch.setattr("time.sleep", lambda s: None)
    return {"parcels": parcels, "catalogue": cat, "tmp": tmp_path}
