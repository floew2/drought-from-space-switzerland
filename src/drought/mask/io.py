"""Writing the purity raster, and site windows in LV95."""

from __future__ import annotations

import numpy as np

from drought.mask.settings import CRS

__all__ = ["site_bounds_lv95", "write_raster"]


def write_raster(bands: list[tuple[np.ndarray, str]], grid: tuple[float, float, float],
                 shape: tuple[int, int], out_tif) -> None:
    """Write purity bands as one float32 GeoTIFF in EPSG:2056.

    Args:
        bands: List of (array, band name) in band order.
        grid: (minx, maxy, res) of the grid.
        shape: (height, width) in cells.
        out_tif: Output path.
    """
    import rasterio
    from rasterio.transform import from_origin

    minx, maxy, res = grid
    height, width = shape
    with rasterio.open(
        out_tif, "w", driver="GTiff", height=height, width=width,
        count=len(bands), dtype="float32", crs=CRS,
        transform=from_origin(minx, maxy, res, res),
        compress="deflate", predictor=3, tiled=True,
    ) as ds:
        for i, (arr, name) in enumerate(bands, 1):
            ds.write(arr.astype(np.float32), i)
            ds.set_band_description(i, name)


def site_bounds_lv95(site: str) -> tuple[float, float, float, float]:
    """A site's bounds in LV95 (EPSG:2056): its AOI polygon's bounding box."""
    from drought.core import aoi

    return aoi.get(site).bounds
