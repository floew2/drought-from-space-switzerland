"""The swissEO backend, offline: small GeoTIFFs stand in for the mosaics and a
fake session for the STAC catalogue."""

from __future__ import annotations

import datetime as dt

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import box

from drought import config, gee, provenance
from drought.config import CFG
from drought.swisseo import extract, grid, scene, stac

X0, Y0, N = 2_720_000.0, 1_230_000.0, 40          # mask: 40 x 40 cells of 10 m
DN = {"b04_10m": 1100, "b08_10m": 3100, "b05_20m": 1600, "b8a_20m": 3100, "b11_20m": 2100}
# reflectance = DN * 1e-4 - 0.1
NDVI = (0.21 - 0.01) / (0.21 + 0.01)
NDMI = (0.21 - 0.11) / (0.21 + 0.11)
NDRE = (0.21 - 0.06) / (0.21 + 0.06)


def _tif(path, arr, x0, y0, res, dtype, desc=None):
    arr = np.asarray(arr)
    arr = arr[None] if arr.ndim == 2 else arr
    with rasterio.open(path, "w", driver="GTiff", width=arr.shape[2], height=arr.shape[1],
                       count=arr.shape[0], dtype=dtype, crs="EPSG:2056",
                       transform=from_origin(x0, y0, res, res)) as ds:
        ds.write(arr.astype(dtype))
        for i, d in enumerate(desc or [], 1):
            ds.set_band_description(i, d)
    return str(path)


def _mosaic(folder, cloudy_rows=0):
    """One mosaic around the mask: 10 m bands on a 10 m grid, 20 m bands on a
    20 m grid whose origin lies OFF the mask's 20 m phase."""
    folder.mkdir(parents=True, exist_ok=True)
    gx, gy, n10 = X0 - 200, Y0 + 200, N + 40
    hrefs = {}
    for k, v in DN.items():
        if k.endswith("10m"):
            hrefs[k] = _tif(folder / f"{k}.tif", np.full((n10, n10), v), gx, gy, 10, "uint16")
        else:
            hrefs[k] = _tif(folder / f"{k}.tif", np.full((n10 // 2 + 1, n10 // 2 + 1), v),
                            gx - 10, gy + 10, 20, "uint16")
    cloud = np.zeros((n10, n10))
    cloud[:cloudy_rows] = 1
    hrefs["cloudmask_10m"] = _tif(folder / "cloud.tif", cloud, gx, gy, 10, "uint8")
    hrefs["terrainmask_10m"] = _tif(folder / "terrain.tif", np.full((n10, n10), 40),
                                    gx, gy, 10, "uint8")
    hrefs["scl_20m"] = _tif(folder / "scl.tif", np.full((n10 // 2 + 1, n10 // 2 + 1), 4),
                            gx - 10, gy + 10, 20, "uint8")
    return hrefs


class _Resp:
    def __init__(self, page): self.page = page
    def raise_for_status(self): pass
    def json(self): return self.page


class _Session:
    """Answers every items query with the mosaics whose date falls in the range."""

    def __init__(self, items): self.items = items

    def get(self, url, params=None, timeout=None):
        a, b = params["datetime"].split("/")
        a, b = a[:10], b[:10]
        feats = [{"id": i, "properties": {"datetime": f"{d}T10:30:00Z"},
                  "assets": {f"x_mosaic_{d}t103000_{k}.tif": {"href": h}
                             for k, h in hrefs.items()}}
                 for i, (d, hrefs) in self.items.items() if a <= d <= b]
        return _Resp({"features": feats, "links": []})


@pytest.fixture
def world(tmp_path, monkeypatch):
    """A tiny site: mask, two mown parcels, one mosaic a year (two in 2026)."""
    from pyproj import Transformer

    monkeypatch.setattr(CFG, "backend", "swisseo")
    for k in ("swisseo", "masks", "labels", "figures"):
        monkeypatch.setattr(CFG.paths, k, str(tmp_path / k))
    monkeypatch.setattr(gee, "_remember", lambda dest, h: None)
    monkeypatch.setattr(CFG.gee.trajectory, "min_pixels", 10)
    monkeypatch.setattr(CFG.gee.cuts, "min_area_ha", 0.5)

    purity = np.full((4, N, N), 0.95, "float32")
    purity[:, :, :5] = 0.5                                     # an impure strip
    _tif(config.mask_path(), purity, X0, Y0, 10, "float32",
         ["mown", "grazed", "intensive", "extensive"])
    # site bbox: the mask, 25 m inside, in WGS84
    bb = Transformer.from_crs(2056, 4326, always_xy=True).transform_bounds(
        X0 + 25, Y0 - N * 10 + 25, X0 + N * 10 - 25, Y0 - 25)
    monkeypatch.setattr(CFG.sites.linthebene, "bbox_wgs84", list(bb))

    parcels = gpd.GeoDataFrame(
        {"parcel_id": [f"{i:020x}" for i in (10, 11)], "lnf_code": [613, 611],
         "area_ha": [1.44, 1.44]},
        geometry=[box(X0 + 60, Y0 - 180, X0 + 180, Y0 - 60),
                  box(X0 + 220, Y0 - 340, X0 + 340, Y0 - 220)], crs=2056)
    items = {}
    for y in extract.years():
        items[f"{y}-07-01"] = (f"{y}-07-01", _mosaic(tmp_path / "img" / str(y)))
    items["2026-07-11"] = ("2026-07-11", _mosaic(tmp_path / "img" / "c", cloudy_rows=40))
    return parcels, _Session(items)


def test_20m_band_is_repeated_and_cut_on_the_10m_window(tmp_path):
    # 20 m cells numbered by column; a window starting 10 m into the first cell
    f = _tif(tmp_path / "c.tif", np.tile(np.arange(10), (10, 1)), 0, 200, 20, "uint16")
    a = scene.read_band(f, grid.Window(10, 190, 6, 2))
    assert a[0].tolist() == [0, 1, 1, 2, 2, 3]


def _stored(it, pl, tmp_site="linthebene"):
    from drought.swisseo import store
    return store.load(store.fetch(it, pl.win, tmp_site), pl.win)


def test_indices_and_masks_from_one_mosaic(world):
    parcels, sess = world
    pl = extract.plan("linthebene", parcels)
    it = stac.search(CFG.sites.linthebene.bbox_wgs84, dt.date(2026, 7, 1),
                     dt.date(2026, 7, 1), session=sess)[0]
    r = extract.process(_stored(it, pl), it.doy, pl)
    ndvi, kndvi, ndmi, ndre, n = r["traj"]["intensive"]
    # stored as int16 x 10 000: equal to 1e-4
    assert ndvi == pytest.approx(NDVI, abs=1e-4)
    assert kndvi == pytest.approx(np.tanh(NDVI ** 2), abs=2e-4)
    assert ndmi == pytest.approx(NDMI, abs=1e-4)
    assert ndre == pytest.approx(NDRE, abs=1e-4)
    assert n == int(pl.regimes["intensive"].sum())        # nothing masked
    # the cloudy mosaic: the top 40 rows of the image are cloud, which covers
    # the top 20 rows of the mask
    cloudy = stac.search(CFG.sites.linthebene.bbox_wgs84, dt.date(2026, 7, 11),
                         dt.date(2026, 7, 11), session=sess)[0]
    r2 = extract.process(_stored(cloudy, pl), cloudy.doy, pl)
    assert 40 < r2["cloud_pct"] < 60
    assert 0 < r2["traj"]["intensive"][4] < n
    assert r2["traj"]["intensive"][2] == pytest.approx(NDMI, abs=1e-4)


def test_stored_file_is_the_whole_window_independent_of_the_mask(world):
    parcels, sess = world
    import rasterio

    from drought.swisseo import store
    pl = extract.plan("linthebene", parcels)
    cloudy = stac.search(CFG.sites.linthebene.bbox_wgs84, dt.date(2026, 7, 11),
                         dt.date(2026, 7, 11), session=sess)[0]
    with rasterio.open(store.fetch(cloudy, pl.win, "linthebene")) as ds:
        a = ds.read(1)
        assert ds.dtypes[0] == "int16" and ds.descriptions == store.BANDS
    # cloud rows empty, everything else kept, grassland or not
    cloud = (a == store.NODATA).all(axis=1)
    assert cloud.any() and not cloud.all()
    assert (a[~cloud] != store.NODATA).all()


def test_run_writes_and_registers_three_products(world):
    parcels, sess = world
    files = extract.run(sites=["linthebene"], session=sess, parcels=parcels, workers=1)
    kinds = {f.name.split("_")[0] for f in files}
    assert {"grassland", "purity", "cuts"} <= kinds
    for f in files:
        kind = ("trajectory" if "trajectory" in f.name else
                "sensitivity" if "sensitivity" in f.name else "cuts")
        rec = provenance.verify(f, kind)
        assert rec["settings"]["backend"] == "swisseo"

    traj = pd.read_csv(next(f for f in files if "trajectory" in f.name))
    assert set(traj.year) == set(extract.years())
    assert (traj.NDMI - NDMI).abs().max() < 1e-4

    from drought.io import load_cuts
    cuts = load_cuts("linthebene", verbose=False)
    assert set(cuts.parcel_id) == set(parcels.parcel_id)
    assert (cuts.NDVI - NDVI).abs().max() < 1e-4

    # a new purity threshold: recomputed from disk, nothing downloaded
    def boom(*a, **k):
        raise AssertionError("downloaded a mosaic again")
    from drought.swisseo import scene as sc
    old = sc.clear_mask
    sc.clear_mask = boom
    try:
        CFG.mask.purity_min, was = 0.4, CFG.mask.purity_min     # takes in the 0.5 strip
        try:
            files2 = extract.run(sites=["linthebene"], session=sess, parcels=parcels, workers=1)
        finally:
            CFG.mask.purity_min = was
    finally:
        sc.clear_mask = old
    t2 = pd.read_csv(next(f for f in files2 if "trajectory" in f.name))
    assert (t2.n_px.to_numpy() > traj.n_px.to_numpy()).all()


def test_earth_engine_hash_and_records_are_unchanged_by_the_backend(monkeypatch):
    monkeypatch.setattr(CFG, "backend", "gee")
    s = gee.render_settings(with_mask=False)
    assert "backend" not in s and "swisseo" not in s
    assert "backend" not in provenance.settings("trajectory")
    h_gee = gee.expected_hash(with_mask=False)
    monkeypatch.setattr(CFG, "backend", "swisseo")
    assert gee.expected_hash(with_mask=False) != h_gee


def test_series_from_another_backend_are_refused(world, monkeypatch):
    parcels, sess = world
    files = extract.run(sites=["linthebene"], session=sess, parcels=parcels, workers=1)
    monkeypatch.setattr(CFG, "backend", "gee")
    with pytest.raises(provenance.ProvenanceError) as exc:
        provenance.verify(files[0], "trajectory")
    assert exc.value.field == "backend"


def test_vhi_on_the_pure_cells_of_each_parcel(world, tmp_path):
    from drought.swisseo import vhi

    parcels, _ = world
    a = np.full((N + 40, N + 40), 30, "uint8")
    a[:20 + 20] = 110                                  # no data over the top 20 mask rows
    f = _tif(tmp_path / "vhi.tif", a, X0 - 200, Y0 + 200, 10, "uint8")
    sess = _Session({"v": ("2026-07-01", {"vegetation-10m": f})})
    out = vhi.run(sites=["linthebene"], parcels=parcels, session=sess, workers=1)
    series = pd.read_csv(out[0])
    per = pd.read_csv(out[1], dtype={"parcel_id": str})
    assert series.vhi_mean.eq(30).all() and (series.n_px > 0).all()
    # the first parcel lies in the no-data rows, the second does not
    assert per.parcel_id.tolist() == [parcels.parcel_id[1]]
    assert per.vhi_mean.tolist() == [30.0]


def test_parcel_id_is_lv95_centroid_and_area():
    from drought.parcel_ids import make_ids
    outer = box(2718000, 1226000, 2718100, 1226200)
    inner = box(2718040, 1226090, 2718060, 1226110)
    g = gpd.GeoSeries([outer.difference(inner), inner], crs=2056)
    # a ring and the field inside it share a centroid; the area tells them apart
    assert make_ids(g) == ["2718050_1226100_19600", "2718050_1226100_400"]
    assert make_ids(g.to_crs(4326)) == make_ids(g)              # any input CRS
