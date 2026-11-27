"""What every product needs for a site: its window, its cells, its parcels.

    c = cells("linthebene")      window, purity per group, the AOI, pure cells
                                 per management system, all grassland cells
    p = parcels("linthebene")    the cadastre (parcels.year) inside the window
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from drought.core import grid
from drought.core.config import CFG, ROOT

__all__ = ["Cells", "cells", "parcels", "rasterize", "write"]

log = logging.getLogger(__name__)


@dataclass
class Cells:
    site: str
    win: grid.Window
    purity: dict[str, np.ndarray]          # group -> purity, 0 outside the mask
    aoi: np.ndarray                         # cells inside the AOI polygon
    regimes: dict[str, np.ndarray]          # system -> pure cells (purity >= purity_min)
    grass: np.ndarray                       # any group purity > 0: the cloud share

    def cloud_pct(self, ix: dict[str, np.ndarray]) -> float:
        """Share (%) of the site's grassland not usable in a mosaic."""
        n = int(self.grass.sum())
        if not n:
            return 100.0
        return float(100.0 * (1 - (self.grass & ~np.isnan(ix["NDVI"])).sum() / n))


def cells(site: str, mask_file=None) -> Cells:
    """The site's cells under the current mask and purity threshold.

    Raises:
        ValueError: If the AOI holds no grassland in the mask (it lies outside it).
    """
    win = grid.site_window(site)
    m = grid.read_mask(win, mask_file)
    inside = grid.aoi_cells(site, win)
    m = {k: np.where(inside, v, 0.0) for k, v in m.items()}
    pmin = float(CFG.mask.purity_min)
    regimes = {r: m[r] >= pmin for r in CFG.analysis.regimes}
    grass = np.zeros(win.shape, bool)
    for a in m.values():
        grass |= a > 0
    if not grass.any():
        raise ValueError(f"{site}: no grassland in the mask here; the mask does not "
                         "cover this AOI (`drought mask` builds it from the parcel data)")
    return Cells(site, win, m, inside, regimes, grass)


def parcels(site: str, win: grid.Window | None = None):
    """The cadastre inside the site window (LV95): parcel_id, lnf_code, area_ha
    (planar, before any edge buffer), geometry."""
    from shapely.geometry import box

    from drought.sources import parcels as cadastre

    g = cadastre.read().to_crs(2056)
    win = win or grid.site_window(site)
    g = g[g.intersects(box(*win.bounds))].reset_index(drop=True)
    g["lnf_code"] = pd.to_numeric(g.lnf_code, errors="coerce")
    g["area_ha"] = g.geometry.area / 1e4
    return g[["parcel_id", "lnf_code", "area_ha", "geometry"]]


def rasterize(geoms, values, win: grid.Window) -> np.ndarray:
    """Label raster: cells whose centre lies inside each polygon get its value."""
    from rasterio.features import rasterize as _r
    shapes = [(g, int(v)) for g, v in zip(geoms, values) if g is not None and not g.is_empty]
    if not shapes:
        return np.zeros(win.shape, np.int32)
    return _r(shapes, out_shape=win.shape, transform=win.transform, fill=0,
              all_touched=False, dtype="int32")


def write(df: pd.DataFrame, path: Path, kind: str) -> Path:
    """Write a product and its provenance record."""
    from drought.core.provenance import register

    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    register(path, kind)
    log.info("  -> %s", path.relative_to(ROOT) if path.is_relative_to(ROOT) else path)
    return path
