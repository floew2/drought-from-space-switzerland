"""Grassland mask: which pixels can be trusted, and for which kind of grassland.

    parcels   reading the annual deliveries; the master parcel set
    purity    the sensor response and PSF-weighted pixel purity
    build     the four-band purity raster, and the area checks
    io        raster writing, site windows
    settings  every mask setting, read once from config.yaml

Figures are in `drought.mask_figures`. The method and its reasoning are in
docs/METHOD.md §2.
"""

from .build import area_table, build_mask, main, output_grid, psf_model_check, verify_mask
from .io import site_bounds_lv95, write_raster
from .parcels import category_history, category_of, load_year, master_parcels, matching_check
from .purity import measured_response_fwhm, psf_caption, psf_params, psf_sigma_m, purity_psf
from .settings import (
                    CRS,
                    GRANDFATHER_CODES,
                    GROUPS,
                    PSF_MODEL,
                    PSF_MTF,
                    PURITY_MIN,
                    RES,
                    SLIVER_FLOOR_M2,
                    SUBPIXEL,
)

__all__ = [
    "CRS", "GRANDFATHER_CODES", "GROUPS", "PSF_MODEL", "PSF_MTF", "PURITY_MIN", "RES",
    "SLIVER_FLOOR_M2", "SUBPIXEL", "area_table", "build_mask", "category_history",
    "category_of", "load_year", "main", "master_parcels", "matching_check",
    "measured_response_fwhm", "output_grid", "psf_caption", "psf_model_check",
    "psf_params", "psf_sigma_m", "purity_psf", "site_bounds_lv95", "verify_mask", "write_raster",
]
