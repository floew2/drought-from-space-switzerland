"""The Earth Engine exports run from Python (drought.engine), tested offline.

Nothing here talks to Earth Engine or Cloud Storage: the server calls are
replaced, and what is tested is the part that decides whether a run may start
and whether a download may enter the analysis.
"""

from __future__ import annotations

import pandas as pd
import pytest

from drought import config, gee, provenance
from drought.engine import cuts, imagery, session, stamp, tasks, trajectory
from drought.provenance import ProvenanceError


@pytest.fixture
def fake_mask(monkeypatch):
    fp = {"path": "data/masks/grassland_mask_sg.tif", "sha256": "a" * 16, "res_m": 10.0,
          "bands": ["mown", "grazed", "intensive", "extensive"], "dtype": "float32",
          "mtime": "2026-10-13T00:00:00"}
    monkeypatch.setattr(provenance, "mask_fingerprint", lambda path=None: dict(fp))
    monkeypatch.setattr(gee, "_remember", lambda dest, h: None)   # no manifest side effect
    return fp


# --------------------------------------------------------------------------- #
# The mask asset must BE the local mask
# --------------------------------------------------------------------------- #
def test_matching_asset_passes(fake_mask):
    assert imagery.check_mask_asset({imagery.MASK_SHA_PROPERTY: "a" * 16}) == "a" * 16


@pytest.mark.parametrize("props,field", [
    ({imagery.MASK_SHA_PROPERTY: "b" * 16}, imagery.MASK_SHA_PROPERTY),   # another build
    ({}, imagery.MASK_SHA_PROPERTY),                                      # uploaded by hand
])
def test_other_or_unstamped_asset_is_refused(fake_mask, props, field):
    with pytest.raises(ProvenanceError) as exc:
        imagery.check_mask_asset(props)
    assert exc.value.field == field


def test_missing_asset_is_refused(fake_mask, monkeypatch):
    monkeypatch.setattr(session, "asset_info", lambda name: None)
    with pytest.raises(ProvenanceError, match="stamp-mask"):
        imagery.check_mask_asset()


def _tiny_mask(path, shift=0.0, swap=False):
    import numpy as np
    import rasterio
    from rasterio.transform import from_origin

    rng = np.random.default_rng(0)
    arr = np.round(rng.uniform(0, 1, (4, 30, 40)) * (rng.uniform(0, 1, (1, 30, 40)) > 0.4),
                   6).astype("float32")
    if swap:
        arr = arr[[1, 0, 2, 3]]
    with rasterio.open(path, "w", driver="GTiff", height=30, width=40, count=4,
                       dtype="float32", crs="EPSG:2056",
                       transform=from_origin(2_722_000 + shift, 1_228_000, 10, 10)) as ds:
        ds.write(arr)
    return path


def test_stamp_compares_the_same_cells(tmp_path):
    """Reading the asset at the sampled cell centres must return the local values."""
    import rasterio

    f = _tiny_mask(tmp_path / "m.tif")
    xy, vals, names = stamp.sample_cells(f, n=60)
    assert names == ["mown", "grazed", "intensive", "extensive"]
    assert (vals.max(axis=1) > 0).sum() >= 45            # three quarters on grassland
    with rasterio.open(f) as ds:                          # the "asset", read at those centres
        remote = {i: list(v) for i, v in enumerate(ds.sample([tuple(p) for p in xy]))}
    assert stamp.compare(vals, remote) == (0, 0.0)


@pytest.mark.parametrize("kind", ["swapped bands", "shifted grid", "cells missing"])
def test_stamp_refuses_a_different_asset(tmp_path, kind):
    import rasterio

    f = _tiny_mask(tmp_path / "m.tif")
    xy, vals, _ = stamp.sample_cells(f, n=60)
    other = _tiny_mask(tmp_path / "o.tif", swap=kind == "swapped bands",
                       shift=10.0 if kind == "shifted grid" else 0.0)
    with rasterio.open(other) as ds:
        remote = {i: list(v) for i, v in enumerate(ds.sample([tuple(p) for p in xy]))}
    if kind == "cells missing":
        remote = {i: v for i, v in remote.items() if i % 7}
    missing, worst = stamp.compare(vals, remote)
    assert missing > 0 or worst > stamp.TOLERANCE


# --------------------------------------------------------------------------- #
# What comes back
# --------------------------------------------------------------------------- #
def test_cut_export_columns_carry_no_geometry():
    sel = cuts.selectors(["20260304_32TNT_NDVI", "20260304_32TNT_NDMI"])
    assert ".geo" not in sel and "system:index" not in sel
    assert sel[:4] == ["parcel_id", "area_ha", "elev_m", "lnf_code"]
    assert sel[4:] == ["20260304_32TNT_NDVI_median", "20260304_32TNT_NDVI_count",
                       "20260304_32TNT_NDMI_median", "20260304_32TNT_NDMI_count"]


def test_cut_years_run_newest_first():
    y = cuts.years()
    assert y[0] == config.CFG.project.current_year and y == sorted(y, reverse=True)


def test_engine_settings_enter_the_config_hash(fake_mask, monkeypatch):
    before = gee.expected_hash()
    monkeypatch.setattr(config.CFG.gee.cuts, "sample_n", 400)
    assert gee.expected_hash() != before


def _trajectory_csv(path):
    years = [*config.CFG.project.baseline_years, config.CFG.project.current_year]
    pd.DataFrame([{"site": "linthebene", "regime": r, "year": y, "date": f"{y}-06-01",
                   "doy": 152, "NDVI": 0.8, "kNDVI": 0.58, "NDMI": 0.35, "NDRE": 0.4,
                   "n_px": 6000} for r in config.CFG.analysis.regimes for y in years]
                 ).to_csv(path, index=False)


def test_fetch_downloads_and_registers(fake_mask, monkeypatch, tmp_path):
    h = tasks.config_hash()
    exp = tasks.Export(None, trajectory.prefix(h), tmp_path, "trajectory")
    monkeypatch.setattr(session, "drive_find", lambda name: [{"id": "x", "name": name}])
    monkeypatch.setattr(session, "drive_download", lambda fid, dest: _trajectory_csv(dest) or dest)
    path = tasks.fetch(exp)
    assert path.name == f"grassland_trajectory_multi_p90_cfg{h}.csv"
    assert provenance.verify(path, "trajectory")["registered_from"] == "earth engine export"


def test_fetch_refuses_an_export_from_other_settings(fake_mask, monkeypatch, tmp_path):
    exp = tasks.Export(None, trajectory.prefix("0" * 12), tmp_path, "trajectory")
    monkeypatch.setattr(session, "drive_find", lambda name: [{"id": "x", "name": name}])
    monkeypatch.setattr(session, "drive_download", lambda fid, dest: _trajectory_csv(dest) or dest)
    with pytest.raises(ProvenanceError):
        tasks.fetch(exp)


# --------------------------------------------------------------------------- #
# One implementation, and machine settings kept out of the method
# --------------------------------------------------------------------------- #
def test_ported_scripts_are_gone():
    for js in ("gee_grassland_trajectory.js", "gee_purity_sensitivity.js",
               "gee_cuts_extract.js"):
        assert not (gee.GEE_DIR / js).exists()
        assert js not in gee.SCRIPTS


def test_local_config_sets_only_where_not_how(tmp_path):
    raw = config.load(local=False).to_dict()
    ok = tmp_path / "config.local.yaml"
    ok.write_text("gee: {project: my-proj, drive_folder: my-folder}\n")
    out = config._overlay_local(raw, ok)
    assert (out["gee"]["project"], out["gee"]["drive_folder"]) == ("my-proj", "my-folder")

    bad = tmp_path / "bad.yaml"
    bad.write_text("mask: {purity_min: 0.5}\n")
    with pytest.raises(ValueError, match="mask.purity_min"):
        config._overlay_local(config.load(local=False).to_dict(), bad)


# --------------------------------------------------------------------------- #
# Every graph-building path runs, against a stand-in for Earth Engine
# --------------------------------------------------------------------------- #
class _Fake:
    """Accepts any attribute or call; `map` calls the function it is given."""

    calls: list = []

    def __init__(self, name="ee"):
        self._name = name

    def __getattr__(self, attr):
        return _Fake(f"{self._name}.{attr}")

    def __call__(self, *args, **kwargs):
        _Fake.calls.append(self._name.rsplit(".", 1)[-1])
        if self._name.endswith(".map"):
            args[0](_Fake("element"))           # run the mapped function once
        if self._name.endswith(".getInfo"):
            return ["20260304_32TNT_NDVI", "20260304_32TNT_NDMI"]
        return _Fake(self._name + "()")

    def __getitem__(self, key):
        return _Fake(f"{self._name}[{key}]")


def test_every_export_builds_against_a_stand_in(fake_mask, monkeypatch):
    from drought.engine import sensitivity

    monkeypatch.setattr(config.CFG.gee, "project", "p")
    monkeypatch.setattr(session, "ee", lambda: _Fake())
    _Fake.calls = []
    traj = trajectory.exports()
    sens = sensitivity.exports()
    cut = cuts.exports()
    assert len(traj) == 1 and traj[0].kind == "trajectory"
    assert len(sens) == 3 * len(config.CFG.analysis.regimes)
    assert len(cut) == len(config.CFG.gee.cuts.sites) * len(cuts.years())
    assert cut[0].selectors[4] == "20260304_32TNT_NDVI_median"
    assert {"reduceRegion", "reduceRegions", "toBands", "linkCollection"} <= set(_Fake.calls)
    smoke = trajectory.exports(smoke=True)[0]
    assert "_SMOKE_cfg" in smoke.prefix and smoke.kind is None


def test_labelled_only_exports_exactly_the_labelled_parcel_years(fake_mask, monkeypatch):
    """The re-score extraction: one export per labelled year, only those parcels."""
    monkeypatch.setattr(config.CFG.gee, "project", "p")
    monkeypatch.setattr(session, "ee", lambda: _Fake())
    labels = {"A|2026": [140], "B|2026": [150], "C|2022": [160]}
    monkeypatch.setattr("drought.labels.load_labels", lambda site: labels)
    assert cuts.labelled("linthebene") == {2026: ["A", "B"], 2022: ["C"]}
    exps = cuts.exports(sites=["linthebene"], labelled_only=True)
    assert [e.prefix.split("_cfg")[0] for e in exps] == [
        "cuts_linthebene_2026_labelled", "cuts_linthebene_2022_labelled"]
    assert all(e.kind == "cuts" and e.dest_dir.name == "labelled" for e in exps)
