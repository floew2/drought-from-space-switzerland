"""The sensor-response model and the purity raster built from it."""

from __future__ import annotations

import numpy as np
import pytest
from shapely import affinity
from shapely.geometry import box

from drought import mask


def test_box_gauss_response_matches_the_stated_derivation():
    p = mask.psf_params(model="box_gauss")
    assert p["mtf_detector"] == pytest.approx(0.637, abs=5e-4)     # sinc(1/2)
    assert p["mtf_optics"] == pytest.approx(0.393, abs=5e-4)       # 0.25 / 0.637
    assert p["sigma_m"] == pytest.approx(4.35, abs=5e-3)
    assert p["fwhm_m"] == pytest.approx(12.7, abs=0.05)
    # The Gaussian alone is narrower than the pixel: never quote it as the response.
    assert 2 * np.sqrt(2 * np.log(2)) * p["sigma_m"] == pytest.approx(10.2, abs=0.05)


def test_single_gaussian_alternative():
    p = mask.psf_params(model="gauss")
    assert p["sigma_m"] == pytest.approx(5.30, abs=5e-3)
    assert p["fwhm_m"] == pytest.approx(12.5, abs=0.05)


def _fields(n=60, seed=0):
    rng = np.random.default_rng(seed)
    return [affinity.rotate(box(x, y, x + rng.uniform(20, 300), y + rng.uniform(15, 60)),
                            rng.uniform(0, 90))
            for x, y in rng.uniform(100, 1900, (n, 2))]


@pytest.mark.parametrize("model", ["box_gauss", "gauss"])
def test_tiled_purity_equals_single_pass(model):
    geoms, grid = _fields(), (0.0, 2000.0, 10.0)
    tiled = mask.purity_psf(geoms, grid, (200, 200), tile=37, model=model)
    whole = mask.purity_psf(geoms, grid, (200, 200), tile=10_000, model=model)
    assert np.array_equal(tiled, whole)


def test_large_field_interior_is_exactly_pure():
    p = mask.purity_psf([box(100, 100, 900, 900)], (0.0, 1000.0, 10.0), (100, 100))
    assert p[50, 50] == 1.0
    assert p.max() <= 1.0 and p.min() >= 0.0


def test_implemented_response_has_the_stated_width():
    """The code must apply box (x) Gaussian, not just the box's sigma.

    Expected widths include the 2.5 m line source used to measure them.
    """
    assert mask.measured_response_fwhm("box_gauss") == pytest.approx(12.76, abs=0.25)
    assert mask.measured_response_fwhm("gauss") == pytest.approx(13.06, abs=0.25)
    # The trap: optical sigma with centre sampling gives a response narrower
    # than either real model, under-weighting the halo.
    assert mask.measured_response_fwhm("gauss", sigma_m=4.352) < 11.3
