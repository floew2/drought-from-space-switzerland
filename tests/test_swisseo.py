"""swissEO end to end, offline: fetch -> store -> products -> runner.

Small GeoTIFFs stand in for the mosaics and a fake catalogue for STAC
(conftest.world).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from conftest import X0, Y0, N, mosaic, ndmi, tif

from drought.core import config, grid, provenance
from drought.core.config import CFG
from drought.sources import fetch, scene, store

SITE = "linthebene"


def _fetch(years=None):
    return fetch.run([SITE], years or config.years(), workers=1)


def _no_download(monkeypatch):
    def boom(*a, **k):
        raise AssertionError("downloaded a mosaic")
    monkeypatch.setattr(scene, "quality", boom)
    monkeypatch.setattr(scene, "read_band", boom)


# --------------------------------------------------------------------------- #
# Reading swissEO
# --------------------------------------------------------------------------- #
def test_20m_band_is_repeated_and_cut_on_the_10m_window(tmp_path):
    # 20 m cells numbered by column; a window starting 10 m into the first cell
    f = tif(tmp_path / "c.tif", np.tile(np.arange(10), (10, 1)), 0, 200, 20, "uint16")
    a = scene.read_band(f, grid.Window(10, 190, 6, 2))
    assert a[0].tolist() == [0, 1, 1, 2, 2, 3]


def test_site_window_lies_on_the_20m_grid_and_ignores_the_mask(world):
    w = grid.site_window(SITE)
    assert w.x0 % 20 == 0 and w.y0 % 20 == 0
    assert w.width * 10 % 20 == 0 and w.height * 10 % 20 == 0


# --------------------------------------------------------------------------- #
# fetch and the store
# --------------------------------------------------------------------------- #
def test_fetch_stores_every_mosaic_and_records_the_catalogue(world):
    assert _fetch() == []
    cur = int(CFG.project.current_year)
    assert len(store.catalogue(SITE, cur)) == 2
    assert all(store.complete(SITE, y) for y in config.years())
    assert store.complete(SITE, cur, "vhi")
    calls = world["catalogue"].calls
    assert _fetch() == []                                   # nothing new to download
    assert world["catalogue"].calls > calls                 # but the catalogue is asked again


def test_dry_run_downloads_nothing(world):
    fetch.run([SITE], config.years(), workers=1, dry_run=True)
    assert not any(store.complete(SITE, y) for y in config.years())


def test_stored_file_keeps_the_quality_layers_and_the_whole_window(world):
    import rasterio
    _fetch([int(CFG.project.current_year)])
    cur = int(CFG.project.current_year)
    with rasterio.open(store.path(SITE, f"{cur}-07-11t103000")) as ds:
        assert tuple(ds.descriptions) == store.BANDS and ds.dtypes[0] == "int16"
        ndmi_band, cloud = ds.read(2), ds.read(4)
    thick = cloud == 1
    assert thick.any() and not thick.all()
    assert (ndmi_band[thick] == store.NODATA).all()          # never usable: not kept
    assert (ndmi_band[~thick] != store.NODATA).all()         # everything else kept


def test_the_cloud_rule_changes_without_a_new_download(world, monkeypatch):
    cur = int(CFG.project.current_year)
    world["catalogue"].add(CFG.swisseo.collection, f"{cur}-07-21t103000", f"{cur}-07-21",
                           mosaic(world["tmp"] / "img" / "thin", cloudy_rows=40,
                                  cloud_class=2))
    _fetch([cur])
    win = grid.site_window(SITE)
    f = store.path(SITE, f"{cur}-07-21t103000")
    strict = store.load(f, win)
    _no_download(monkeypatch)
    monkeypatch.setattr(CFG.swisseo.usable, "cloud_clear", [0, 2])
    lenient = store.load(f, win)
    assert np.isnan(strict["NDMI"]).sum() > 0
    assert np.isnan(lenient["NDMI"]).sum() == 0


def test_a_read_cut_off_in_transit_is_retried(world, monkeypatch):
    from rasterio.errors import RasterioIOError
    real, calls = scene.quality, []

    def flaky(item, win):
        calls.append(item.id)
        if len(calls) == 1:
            raise RasterioIOError("TIFFReadEncodedTile() failed")
        return real(item, win)

    monkeypatch.setattr(scene, "quality", flaky)
    assert _fetch([int(CFG.project.current_year)]) == []


def test_a_mosaic_that_keeps_failing_is_reported_and_blocks_the_products(world, monkeypatch):
    from rasterio.errors import RasterioIOError

    from drought.products import condition
    real = scene.quality

    def broken_2020(item, win):
        if item.id.startswith("2020"):
            raise RasterioIOError("TIFFReadEncodedTile() failed")
        return real(item, win)

    monkeypatch.setattr(scene, "quality", broken_2020)
    assert _fetch() == ["2020-07-01t103000"]
    with pytest.raises(FileNotFoundError, match="drought fetch"):
        condition.run(SITE)
    monkeypatch.setattr(scene, "quality", real)             # the network is back
    assert _fetch() == []
    assert condition.run(SITE).exists()


# --------------------------------------------------------------------------- #
# Products
# --------------------------------------------------------------------------- #
def test_condition_series(world):
    from drought.products import condition
    _fetch()
    df = pd.read_csv(condition.run(SITE))
    assert set(df.year) == set(config.years())
    cur = int(CFG.project.current_year)
    base = df[(df.year == 2019) & (df.regime == "intensive")].NDMI.iloc[0]
    assert base == pytest.approx(ndmi(3100, 2100), abs=1e-4)
    now = df[(df.year == cur) & (df.regime == "intensive")].sort_values("doy")
    assert (now.NDMI < base).all()
    assert now.n_px.iloc[1] < now.n_px.iloc[0]               # the cloudy date has fewer cells
    assert provenance.verify(condition.output(SITE), "condition")


def test_cut_series_sample_and_labelled(world):
    import json

    from drought.analysis.io import load_cuts
    from drought.products import cuts
    cur = int(CFG.project.current_year)
    pid = world["parcels"].parcel_id[0]
    config.labels_path(SITE).write_text(json.dumps({f"{pid}|{cur}": [180]}))
    _fetch()
    files = cuts.run(SITE)
    assert cuts.output(SITE, cur, labelled=True) in files
    df = load_cuts(SITE, verbose=False)
    assert set(df.parcel_id) == set(world["parcels"].parcel_id)
    lab = load_cuts(SITE, verbose=False, labelled=True)
    assert set(lab.parcel_id) == {pid}


def test_sensitivity_anomaly_and_vhi(world):
    from drought.products import anomaly, sensitivity, vhi
    _fetch()
    for f in sensitivity.run(SITE):
        assert provenance.is_current(f, "sensitivity")
    a = pd.read_csv(anomaly.run(SITE))
    assert (a.z.dropna() < -1).all()                         # the current year is dry
    series, parcels = vhi.run(SITE)
    s = pd.read_csv(series)
    assert s.vhi_mean.eq(25).all()
    p = pd.read_csv(parcels, dtype={"parcel_id": str})
    # the first parcel lies in the no-data rows, the second does not
    assert p.parcel_id.tolist() == [world["parcels"].parcel_id[1]]


def test_a_new_purity_threshold_recomputes_without_downloading(world, monkeypatch):
    from drought.products import condition
    _fetch()
    before = pd.read_csv(condition.run(SITE))
    _no_download(monkeypatch)
    monkeypatch.setattr(CFG.mask, "purity_min", 0.4)        # takes in the 0.5 strip
    assert not provenance.is_current(condition.output(SITE), "condition")
    after = pd.read_csv(condition.run(SITE))
    assert (after.n_px.to_numpy() > before.n_px.to_numpy()).all()


def test_geojson_aoi_is_a_site_and_its_polygon_limits_the_cells(world, tmp_path):
    import geopandas as gpd
    from shapely.geometry import box

    from drought.core import aoi
    from drought.products.common import cells
    f = tmp_path / "half.geojson"
    gpd.GeoDataFrame(geometry=[box(X0, Y0 - N * 10, X0 + N * 5, Y0)], crs=2056) \
        .to_crs(4326).to_file(f, driver="GeoJSON")
    name = aoi.resolve(str(f))
    assert name == "half" and CFG.sites.half.label == "half"
    c = cells(name)
    full = cells(SITE)
    assert 0 < c.regimes["intensive"].sum() < full.regimes["intensive"].sum()


# --------------------------------------------------------------------------- #
# The runner
# --------------------------------------------------------------------------- #
def test_runner_runs_what_is_needed_once(world):
    from drought.pipeline import runner
    plan = [(s, w) for s, _site, w in runner.plan("condition", [SITE])]
    assert plan == [("fetch", True), ("mask", False), ("condition", True)]
    ran = runner.run("condition", [SITE], workers=1)
    assert ran == [("fetch", SITE), ("condition", SITE)]
    assert runner.run("condition", [SITE]) == []             # everything up to date
    assert runner.run("condition", [SITE], force=True) == [("condition", SITE)]


def test_runner_reruns_a_stage_whose_setting_changed(world, monkeypatch):
    from drought.pipeline import runner
    runner.run("condition", [SITE], workers=1)
    monkeypatch.setattr(CFG.swisseo.usable, "cloud_clear", [0, 2])
    assert runner.run("condition", [SITE]) == [("condition", SITE)]


def test_window_cells_match_the_mask_grid(world):
    """The window starts on the 20 m grid, the mask on its own 10 m grid: their cells coincide."""
    w = grid.site_window(SITE)
    m = grid.read_mask(w)
    assert m["intensive"].shape == w.shape
    col0 = int(round((X0 - w.x0) / 10))
    row0 = int(round((w.y0 - Y0) / 10))
    assert m["intensive"][row0 + 1, col0 + 10] == pytest.approx(0.95)
    assert m["intensive"][row0 + 1, col0 + 2] == pytest.approx(0.5)
