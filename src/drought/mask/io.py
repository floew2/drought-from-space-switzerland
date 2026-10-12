"""Writing the purity raster, and site windows in LV95."""

from __future__ import annotations

import numpy as np

from .settings import CRS

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
    """A site's bounding box from config, transformed WGS84 -> LV95.

    Args:
        site: Site key in config.yaml.

    Returns:
        (x0, y0, x1, y1) in EPSG:2056.
    """
    from pyproj import Transformer

    from ..config import site_cfg

    b = site_cfg(site).bbox_wgs84
    t = Transformer.from_crs("EPSG:4326", CRS, always_xy=True)
    x0, y0 = t.transform(b[0], b[1])
    x1, y1 = t.transform(b[2], b[3])
    return min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1)
