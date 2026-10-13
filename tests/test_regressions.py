"""One test per bug that actually happened (METHOD §6).

Where it can, each test also runs the OLD logic and shows that it fails the
same check. That is the evidence the check would have caught the bug, rather
than merely passing on the fixed code.
"""

from __future__ import annotations

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import yaml
from shapely import affinity
from shapely.geometry import Polygon, box

from drought import config, detect, mask
from drought.mask.parcels import _overlap_years, _repair
from drought.mask.purity import _burn

CRS = "EPSG:2056"
X0, Y0 = 2_722_000.0, 1_227_500.0


def _fields(n=40, seed=1):
    """Long, narrow fields, rotated and partly abutting: Swiss grassland in miniature."""
    rng = np.random.default_rng(seed)
    out, x = [], X0
    for _ in range(n):
        w, h = rng.uniform(80, 250), rng.uniform(15, 45)
        out.append(affinity.rotate(box(x, Y0, x + w, Y0 + h), 17, origin=(X0, Y0)))
        x += w + (0 if rng.uniform() < 0.6 else rng.uniform(3, 12))
    return out


# --------------------------------------------------------------------------- #
# Bug 1: persistence applied per cell zeroed the halo                (METHOD §6)
# --------------------------------------------------------------------------- #
def test_area_at_purity_zero_exceeds_polygon_area():
    """The response reaches past every boundary, so the ratio at > 0 must exceed 1.

    The old mask multiplied purity by a per-cell "grassland every year" layer.
    That zeroed the halo outside the fields and pushed the ratio below 1,
    which is how the bug was found.
    """
    geoms = _fields()
    b = gpd.GeoSeries(geoms).total_bounds
    grid = (np.floor(b[0] / 10) * 10 - 50, np.ceil(b[3] / 10) * 10 + 50, 10.0)
    shape = (int((grid[1] - b[1]) / 10) + 10, int((b[2] - grid[0]) / 10) + 10)
    polygon_ha = sum(g.area for g in geoms) / 1e4

    purity = mask.purity_psf(geoms, grid, shape)
    assert (purity > 0).sum() * 0.01 / polygon_ha > 1.0

    # The old design: a binary per-cell persistence layer multiplied into the raster.
    persistent_cells = _burn(geoms, grid, shape)
    old = purity * persistent_cells
    assert (old > 0).sum() * 0.01 / polygon_ha < 1.0


# --------------------------------------------------------------------------- #
# Bug 3: overlap matching silently dropped rows for a subset          (METHOD §6)
# --------------------------------------------------------------------------- #
def _frames():
    geoms = _fields(20)
    cur = gpd.GeoDataFrame({"_code": [613] * 20}, geometry=geoms, crs=CRS)
    # The past delivery: the same ground, redrawn slightly larger.
    past = gpd.GeoDataFrame({"_code": [613] * 20},
                            geometry=[g.buffer(0.5) for g in geoms], crs=CRS)
    return cur, past


def test_overlap_matching_on_a_subset_keeps_every_row():
    """A per-category call passes a SUBSET of rows; every one must still match.

    The bug: a frame built from {"_i": cur.index} next to a separately-indexed
    GeoSeries. pandas gave the column a fresh RangeIndex, the geometry kept its
    own, and a subset lost rows silently: 36 % matched where 100 % was right,
    while the full set passed.
    """
    cur, past = _frames()
    subset = cur.iloc[1::2]                      # non-contiguous index: 1, 3, 5, ...
    frac = _overlap_years(subset, past, [613], 0.8)
    assert list(frac.index) == list(subset.index)
    assert (frac >= 0.99).all()

    # The old construction, on the same subset: half the rows lose their
    # geometry and drop out of the join, with no error.
    old = gpd.GeoDataFrame({"_i": subset.index}, geometry=subset.geometry, crs=CRS)
    pairs = gpd.sjoin(old[old.geometry.notna()], past[["geometry"]], predicate="intersects")
    assert pairs["_i"].nunique() < len(subset)


# --------------------------------------------------------------------------- #
# Geometry repair: pathological rings must survive, and keep their area
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("name,ring,area", [
    ("self-touching ring", [(0, 0), (4, 0), (4, 4), (2, 2), (2, 4), (0, 4), (2, 2), (0, 0)], 10.0),
    ("figure-eight (bowtie)", [(0, 0), (4, 4), (4, 0), (0, 4), (0, 0)], 8.0),
    ("ring touching its own edge", [(0, 0), (10, 0), (10, 10), (5, 0), (0, 10), (0, 0)], 50.0),
])
def test_geometry_repair_survives_pathological_rings(name, ring, area):
    """A self-touching ring once threw a side-location conflict minutes into a run.

    Repair must make it valid, keep all of its area, and let the overlay run.
    buffer(0), which the repair used to try first, drops one lobe of a
    figure-eight without complaint.
    """
    poly = affinity.translate(Polygon(ring), X0, Y0)
    assert not poly.is_valid
    g = _repair(gpd.GeoDataFrame({"_code": [613]}, geometry=[poly], crs=CRS))
    assert len(g) == 1 and g.geometry.is_valid.all()
    assert g.geometry.iloc[0].area == pytest.approx(area)
    assert g.geom_type.isin(["Polygon", "MultiPolygon"]).all()
    _overlap_years(g, g, [613], 0.8)                     # the overlay must not raise

    if name.startswith("figure-eight"):
        assert poly.buffer(0).area == pytest.approx(area / 2)   # the old first step


# --------------------------------------------------------------------------- #
# The detector finds planted cuts, and only those
# --------------------------------------------------------------------------- #
def _series(cuts=(), decline=None, dip=None, step=5):
    """One parcel-year: green-up, plateau, cuts with a week of stubble and ~3 weeks regrowth."""
    days = np.arange(60, 306, step).astype(float)
    env = 0.06 + 0.36 * np.clip((days - 75) / 50, 0, 1)
    if decline:                                  # (start, end, total drop)
        a, b, d = decline
        env -= d * np.clip((days - a) / (b - a), 0, 1)
    ndmi = env.copy()
    for c in cuts:
        dt = days - c
        removed = 0.32 * np.clip(1 - (dt - 6) / 22, 0, 1)
        removed[dt < 0] = 0
        ndmi = np.minimum(ndmi, env - removed)
    if dip is not None:                          # one bad date that bounces straight back
        ndmi[np.argmin(abs(days - dip))] -= 0.28
    ndvi = 0.90 - 0.45 * np.exp(-5 * ndmi)
    if dip is not None:
        ndvi[np.argmin(abs(days - dip))] -= 0.2
    return pd.DataFrame({"doy": days, "NDMI": ndmi, "NDVI": ndvi, "lnf_code": 613})


def test_detector_finds_planted_cuts():
    planted = [130, 170, 212, 255]
    found = detect.detect_parcel_year(_series(planted))["cut_dates"]
    assert len(found) == len(planted)
    # dated at the first observation clearly down: within a week after the cut
    assert all(0 <= f - p <= 7 for f, p in zip(sorted(found), planted))


def test_detector_ignores_a_drought_decline():
    """Mowing is a step, drought a slope. Without the drop-rate rule the detector
    walked down a decline and called each stage a cut.

    0.30 over 35 days (0.0086 / day) is the case only that rule can reject: it
    exceeds the 0.20 minimum drop within the 26-day lookback, but falls more
    slowly than 0.010 / day.
    """
    s = _series(decline=(185, 220, 0.30))
    assert detect.detect_parcel_year(s)["cut_dates"] == []
    assert detect.detect_parcel_year(s, min_drop_rate=0)["cut_dates"], "rule off: fires"


def test_detector_ignores_a_single_date_dip():
    """Cut grass stays short for weeks; a dip that bounces straight back is atmosphere."""
    s = _series(dip=160)
    assert detect.detect_parcel_year(s)["cut_dates"] == []
    assert detect.detect_parcel_year(s, persist_frac=-1e9)["cut_dates"], "rule off: fires"


# --------------------------------------------------------------------------- #
# Config round-trip
# --------------------------------------------------------------------------- #
def test_config_round_trip():
    """config.yaml -> Config -> dict is lossless, and switching profiles restores it exactly.

    The settings travel as one object, and provenance dumps that object
    beside every result. If the round trip lost or altered a value, the
    provenance would describe a configuration that never ran.
    """
    raw = yaml.safe_load((config.ROOT / "config.yaml").read_text())
    assert config.load().to_dict() == raw
    assert yaml.safe_load(yaml.safe_dump(config.CFG.to_dict())) == config.CFG.to_dict()
    assert config.CFG.get("detect.min_drop_abs") == raw["detect"]["min_drop_abs"]
    assert config.CFG.get("no.such.key", "absent") == "absent"

    before = config.CFG.to_dict()
    try:
        config.activate("synthetic")
        assert config.CFG.to_dict() != before
        config.activate("real")
        assert config.CFG.to_dict() == before
    finally:
        config.activate("real")
    assert config.provenance()["config"] == before


# --------------------------------------------------------------------------- #
# Bug 13: the build summary divided clipped pixels by unclipped polygons
# --------------------------------------------------------------------------- #
def test_area_ratio_reference_is_clipped_to_the_raster():
    """With mask.clip_bbox_lv95 set, parcels reach past the grid.

    The first real build printed ratios at > 0 of 0.74 (grazed): summering
    pastures lie mostly outside the clip box, so their polygon area counted
    in the denominator while their pixels could not exist. Only the build's
    own summary was affected; area_table always clipped.
    """
    from drought.mask.build import window_polygon_ha

    geoms = _fields()
    # a large field that lies mostly beyond the grid's right-hand edge
    far = box(X0 + 4_000, Y0 + 200, X0 + 9_000, Y0 + 400)
    geoms = [*geoms, far]
    grid = (X0 - 100, Y0 + 600, 10.0)
    shape = (90, int((X0 + 4_500 - grid[0]) / 10))          # ends 500 m into `far`

    purity = mask.purity_psf(geoms, grid, shape)
    observed = (purity > 0).sum() * 0.01
    assert observed / window_polygon_ha(geoms, grid, shape) > 1.0

    unclipped = sum(g.area for g in geoms) / 1e4             # the old denominator
    assert observed / unclipped < 1.0
