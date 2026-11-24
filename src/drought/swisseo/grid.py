"""The site windows, on the mask's 10 m LV95 grid, and the mask inside them.

swissEO's 10 m bands lie on a grid whose lines fall on multiples of 10 m in
LV95, as do the mask's. A site window is therefore cut from the mask grid and
read from every 10 m band without resampling; the 20 m bands are read on their
own grid and each cell is repeated 2 x 2 (nearest neighbour, as Earth Engine
does at scale 10).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..config import CFG, mask_path

__all__ = ["Window", "site_window", "read_mask"]


@dataclass(frozen=True)
class Window:
    """A rectangle of whole 10 m cells: upper-left corner and shape."""
    x0: float
    y0: float               # TOP edge (LV95 northing)
    width: int
    height: int
    res: float = 10.0

    @property
    def bounds(self) -> tuple[float, float, float, float]:
        return (self.x0, self.y0 - self.height * self.res,
                self.x0 + self.width * self.res, self.y0)

    @property
    def shape(self) -> tuple[int, int]:
        return (self.height, self.width)

    @property
    def transform(self):
        from rasterio.transform import from_origin
        return from_origin(self.x0, self.y0, self.res, self.res)


def _mask_meta(path=None):
    import rasterio
    with rasterio.open(path or mask_path()) as ds:
        return ds.transform, ds.width, ds.height, list(ds.descriptions)


def site_window(site: str, path=None) -> Window:
    """The site rectangle (config bbox_wgs84) in LV95, snapped onto the mask grid
    and clipped to the mask's extent."""
    from pyproj import Transformer

    tr, w, h, _ = _mask_meta(path)
    res = tr.a
    bb = getattr(CFG.sites, site).bbox_wgs84
    x0, y0, x1, y1 = Transformer.from_crs(4326, 2056, always_xy=True).transform_bounds(*bb)
    c0 = max(int(np.floor((x0 - tr.c) / res)), 0)
    c1 = min(int(np.ceil((x1 - tr.c) / res)), w)
    r0 = max(int(np.floor((tr.f - y1) / res)), 0)
    r1 = min(int(np.ceil((tr.f - y0) / res)), h)
    if c1 <= c0 or r1 <= r0:
        raise ValueError(f"site {site!r} does not overlap the mask")
    return Window(tr.c + c0 * res, tr.f - r0 * res, c1 - c0, r1 - r0, res)


def read_mask(win: Window, path=None) -> dict[str, np.ndarray]:
    """Purity per management group inside the window, by band name."""
    import rasterio
    from rasterio.windows import from_bounds

    with rasterio.open(path or mask_path()) as ds:
        from ..mask.settings import GROUPS
        names = [d or g for d, g in zip(ds.descriptions, GROUPS)]
        rw = from_bounds(*win.bounds, transform=ds.transform).round_offsets().round_lengths()
        arr = ds.read(window=rw, boundless=True, fill_value=0).astype("float32")
    if arr.shape[1:] != win.shape:
        raise ValueError(f"mask window {arr.shape[1:]} != {win.shape}")
    return {n: np.nan_to_num(arr[i]) for i, n in enumerate(names)}
