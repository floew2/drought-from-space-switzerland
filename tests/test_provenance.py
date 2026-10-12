"""Stale files are refused, and the refusal names what is stale (METHOD §6, bug 6).

The worked example from the project's history: a trajectory export with the
right schema but the wrong vintage, carrying management categories (mown /
grazed) that the pipeline no longer uses, and a mask from an older build.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from drought import config, gee, provenance
from drought.provenance import ProvenanceError, record, register, verify, write_sidecar


def _trajectory_csv(path, regimes=("intensive", "extensive")):
    years = [*config.CFG.project.baseline_years, config.CFG.project.current_year]
    rows = [{"site": "linthebene", "regime": r, "year": y, "date": f"{y}-06-01", "doy": 152,
             "NDVI": 0.8, "kNDVI": 0.58, "NDMI": 0.35, "NDRE": 0.4, "n_px": 6000}
            for r in regimes for y in years]
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


@pytest.fixture
def fake_mask(monkeypatch):
    """A current mask fingerprint without building a mask."""
    fp = {"path": "data/masks/grassland_mask_sg.tif", "sha256": "a" * 16, "res_m": 10.0,
          "bands": ["mown", "grazed", "intensive", "extensive"], "dtype": "float32",
          "mtime": "2026-10-01T00:00:00"}
    monkeypatch.setattr(provenance, "mask_fingerprint", lambda path=None: dict(fp))
    return fp


def test_current_file_passes(tmp_path, fake_mask):
    f = _trajectory_csv(tmp_path / "traj.csv")
    register(f, "trajectory")
    assert verify(f, "trajectory")["content"]["regimes"] == ["extensive", "intensive"]


def test_missing_record_is_refused(tmp_path, fake_mask):
    f = _trajectory_csv(tmp_path / "traj.csv")
    with pytest.raises(ProvenanceError, match="provenance record"):
        verify(f, "trajectory")


def test_old_categories_are_refused_at_registration(tmp_path, fake_mask):
    """The CSV that was nearly published: mown / grazed, not intensive / extensive."""
    f = _trajectory_csv(tmp_path / "grassland_trajectory_multi.csv", ("mown", "grazed"))
    with pytest.raises(ProvenanceError) as exc:
        register(f, "trajectory")
    assert exc.value.field == "content.regimes"
    assert exc.value.recorded == ["grazed", "mown"]


def test_record_from_other_settings_is_refused_naming_the_setting(tmp_path, fake_mask):
    f = _trajectory_csv(tmp_path / "traj.csv")
    stale = record("trajectory", f, provenance.summarise("trajectory", f))
    stale["settings"]["mask.purity_min"] = 0.8               # produced at another threshold
    write_sidecar(f, stale)
    with pytest.raises(ProvenanceError) as exc:
        verify(f, "trajectory")
    assert exc.value.field == "mask.purity_min"
    assert (exc.value.recorded, exc.value.current) == (0.8, 0.9)


def test_record_from_another_mask_is_refused(tmp_path, fake_mask):
    f = _trajectory_csv(tmp_path / "traj.csv")
    register(f, "trajectory")
    side = provenance.sidecar_path(f)
    rec = json.loads(side.read_text())
    rec["mask"]["sha256"] = "b" * 16                         # an older mask build
    side.write_text(json.dumps(rec))
    with pytest.raises(ProvenanceError) as exc:
        verify(f, "trajectory")
    assert exc.value.field == "mask.sha256"


def test_file_changed_after_recording_is_refused(tmp_path, fake_mask):
    f = _trajectory_csv(tmp_path / "traj.csv")
    register(f, "trajectory")
    with open(f, "a") as fh:
        fh.write("linthebene,intensive,2026,2026-09-30,273,0.8,0.58,0.4,0.4,6000\n")
    with pytest.raises(ProvenanceError) as exc:
        verify(f, "trajectory")
    assert exc.value.field == "file_sha256"


def test_a_20m_uint8_mask_is_not_a_current_mask(tmp_path, monkeypatch):
    """The mask that sat beside the stale CSV: 20 m, uint8. Not what the pipeline builds."""
    import rasterio
    from rasterio.transform import from_origin

    f = tmp_path / "grassland_mask_sg.tif"
    with rasterio.open(f, "w", driver="GTiff", height=4, width=4, count=1, dtype="uint8",
                       crs="EPSG:2056", transform=from_origin(2_700_000, 1_250_000, 20, 20)) as ds:
        ds.write(np.ones((4, 4), np.uint8), 1)
    write_sidecar(f, record("mask", f))
    with pytest.raises(ProvenanceError, match="resolution / dtype"):
        verify(f, "mask")


def test_require_false_warns_instead(tmp_path, fake_mask, monkeypatch):
    f = _trajectory_csv(tmp_path / "traj.csv")
    monkeypatch.setattr(config.CFG.provenance, "require", False)
    with pytest.warns(UserWarning, match="provenance not verified"):
        verify(f, "trajectory")


def test_export_named_with_another_config_hash_is_refused(tmp_path, fake_mask):
    f = _trajectory_csv(tmp_path / "grassland_trajectory_multi_p90_cfg000000000000.csv")
    with pytest.raises(ProvenanceError, match="Re-render"):
        register(f, "trajectory")


# --------------------------------------------------------------------------- #
# The Earth Engine scripts
# --------------------------------------------------------------------------- #
def test_committed_gee_scripts_match_the_renderer():
    """CI runs this: a hand edit to a generated value, or a config change not
    rendered, fails the build."""
    assert gee.check() == []


def test_committed_scripts_carry_no_project_id():
    """Every Earth Engine project named in a committed script is the placeholder,
    and no asset path is written out by hand."""
    import re

    for script in gee.SCRIPTS:
        text = (gee.GEE_DIR / script).read_text()
        assert set(re.findall(r"var PROJECT = '([^']*)'", text.replace('"', "'"))) <= {
            gee.PLACEHOLDER}
        assert re.findall(r"projects/[A-Za-z0-9_-]+/assets", text) == []
        assert gee.PLACEHOLDER in text or "ASSETS" not in gee.SCRIPTS[script]


def test_paste_ready_render_needs_a_project(tmp_path):
    with pytest.raises(ValueError, match="gee.project"):
        gee.render(out=tmp_path)


# --------------------------------------------------------------------------- #
# results.json never mixes vintages
# --------------------------------------------------------------------------- #
def test_results_drop_and_refuse_stale_sections(tmp_path, monkeypatch, fake_mask):
    from drought import results

    monkeypatch.setattr(results, "path", lambda: tmp_path / "results.json")
    results.update("condition", {"linthebene/intensive": {"integral_pct": -33.0}})
    assert "condition" in results.load()["sections"]

    # The mask is rebuilt: the earlier section no longer describes the current state.
    fake_mask["sha256"] = "c" * 16
    with pytest.raises(RuntimeError, match="condition"):
        results.load()
    results.update("precipitation", {"precip_mm": 501.0})       # writing drops the stale one
    assert set(results.load()["sections"]) == {"precipitation"}
