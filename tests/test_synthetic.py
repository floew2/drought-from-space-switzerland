"""The synthetic data is calibrated, deterministic and kept apart from real data.

These cover the generator and the plumbing it depends on. The five tests tied
to bugs that occurred in the real analysis live in their own modules.
"""

from __future__ import annotations

import numpy as np
import pytest

from drought import config
from drought.synthetic import parcels as P
from drought.synthetic import targets as T

SCALE = 0.125


@pytest.fixture(scope="module")
def population():
    return P.draw_population(np.random.default_rng(7), SCALE)


def _shape_index(w, h):
    return 2 * (w + h) / np.sqrt(w * h)


def test_aspect_for_shape_inverts_the_shape_index():
    s = np.array([4.0, 4.5, 5.72, 6.29, 10.0, 16.0])
    r = P.aspect_for_shape(s)
    assert np.allclose(_shape_index(r, 1.0), s)
    # The formula an earlier attempt used, r = s^2 / 16, lands near 4.2 for a
    # target of 5.7 — the reason the derivation is spelled out in the module.
    assert _shape_index(5.72**2 / 16, 1.0) < 4.3


def test_area_fit_hits_median_small_share_and_mean():
    tgt = T.GEOMETRY["intensive"]
    n, ha = T.CATEGORY_TOTALS["intensive"]
    mu, sigma, b = P.fit_area(tgt["median_ha"], tgt["lt_05"], ha / n)
    assert P._tln_cdf(tgt["median_ha"], mu, sigma, b) == pytest.approx(0.5, abs=1e-6)
    assert P._tln_cdf(0.5, mu, sigma, b) == pytest.approx(tgt["lt_05"], abs=1e-6)
    assert P._tln_mean(mu, sigma, b) == pytest.approx(ha / n, rel=1e-6)


@pytest.mark.parametrize("cat,components", [
    ("intensive", {"I"}), ("extensive", {"Em", "Eg"}), ("grazed", {"Eg", "Go", "S"}),
])
def test_population_matches_the_real_geometry(population, cat, components):
    # The targets describe the MASTER set: the parcels planted to survive it.
    master = population[population.status.isin(["stable", "exempt_930"])]
    sub = master[master.component.isin(components)]
    area = sub.w * sub.h / 1e4
    shape = _shape_index(sub.w, sub.h)
    tgt = T.GEOMETRY[cat]
    assert np.median(area) == pytest.approx(tgt["median_ha"], rel=0.05)
    assert (area < 0.5).mean() == pytest.approx(tgt["lt_05"], abs=0.03)
    assert np.median(shape) == pytest.approx(tgt["median_shape"], abs=0.1)
    assert area.min() >= P.MIN_HA


def test_shape_index_has_the_long_upper_tail(population):
    shape = _shape_index(population.w, population.h)
    assert 0.03 < (shape > 10).mean() < 0.20           # a minority exceed 10
    # the one real tail statistic: 2.3 % exceed 16 (narrow strips)
    assert (shape > 16).mean() == pytest.approx(T.SHAPE_TAIL_GT16, abs=0.005)


def test_generation_is_deterministic():
    a, da = P.generate(np.random.default_rng(3), 0.02)
    b, db = P.generate(np.random.default_rng(3), 0.02)
    assert a.equals(b)
    assert all(da[y].geometry.geom_equals_exact(db[y].geometry, 0).all() for y in da)


def test_deliveries_redraw_ids_and_hold_930_back(population):
    pop, dl = P.generate(np.random.default_rng(5), 0.02)
    assert set(dl) == set(P.YEARS)
    # object IDs never carry over between deliveries
    assert not set(dl[2024].objectid) & set(dl[2025].objectid)
    # Sömmerung exists only in the latest delivery
    assert (dl[2025].lnf_code == 930).any()
    assert not any((dl[y].lnf_code == 930).any() for y in P.YEARS[:-1])
    assert all(g.crs.to_epsg() == 2056 for g in dl.values())
    assert all(g.geometry.is_valid.all() for g in dl.values())


def test_synthetic_profile_never_shares_a_folder_with_real_data():
    try:
        config.activate("synthetic")
        syn = {k: v for k, v in vars(config.CFG.paths).items()}
        assert all(v.startswith(("data/synthetic/", "output/synthetic/")) for v in syn.values())
        config.activate("real")
        real = vars(config.CFG.paths)
        assert not any("synthetic" in v for v in real.values())
        assert set(syn) == set(real), "every real path needs a synthetic counterpart"
    finally:
        config.activate("real")


def test_every_saved_figure_is_stamped(tmp_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from drought.synthetic import STAMP, stamp_figures

    stamp_figures()
    fig, _ = plt.subplots()
    fig.savefig(tmp_path / "x.png")
    assert any(t.get_text() == STAMP for t in fig.texts)
    plt.close(fig)


def test_station_series_is_calibrated():
    from drought.synthetic import station

    df = station.generate(np.random.default_rng(11))
    end = T.PRECIP["end_doy"]
    p26 = df.rre150d0[(df.index.year == 2026) & (df.index.dayofyear <= end)].sum()
    assert p26 == pytest.approx(T.PRECIP["p_2026_mm"], abs=1.5)
    assert df.index.max().dayofyear == end
