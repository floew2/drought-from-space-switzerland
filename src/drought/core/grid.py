"""Site windows on the 10 m LV95 grid, and the mask and AOI inside them.

swissEO's 10 m bands lie on a grid whose lines fall on multiples of 10 m in
LV95, as do the mask's. A site window is the AOI's bounding box widened to
multiples of 20 m, so it lines up with the 10 m and the 20 m bands alike and
does not depend on the mask: the same window, and therefore the same stored
mosaics, serve any mask. The 10 m bands are read without resampling; the 20 m
bands on their own grid, each cell repeated 2 x 2 (nearest neighbour).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from drought.core.config import mask_path

__all__ = ["Window", "aoi_cells", "read_mask", "site_window"]


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


SNAP_M = 20.0


def site_window(site: str) -> Window:
    """The AOI's bounding box in LV95, widened to multiples of SNAP_M."""
    from drought.core import aoi

    x0, y0, x1, y1 = aoi.get(site).bounds
    x0, y0 = np.floor(x0 / SNAP_M) * SNAP_M, np.floor(y0 / SNAP_M) * SNAP_M
    x1, y1 = np.ceil(x1 / SNAP_M) * SNAP_M, np.ceil(y1 / SNAP_M) * SNAP_M
    res = 10.0
    return Window(float(x0), float(y1), int(round((x1 - x0) / res)),
                  int(round((y1 - y0) / res)), res)


def aoi_cells(site: str, win: Window) -> np.ndarray:
    """Cells whose centre lies inside the AOI polygon (all, for a rectangle)."""
    from rasterio.features import rasterize

    from drought.core import aoi

    return rasterize([(aoi.get(site).geometry, 1)], out_shape=win.shape,
                     transform=win.transform, fill=0, dtype="uint8").astype(bool)


def read_mask(win: Window, path=None) -> dict[str, np.ndarray]:
    """Purity per management group inside the window, by band name; 0 outside
    the mask's extent."""
    import rasterio
    from rasterio.windows import from_bounds

    with rasterio.open(path or mask_path()) as ds:
        from drought.mask.settings import GROUPS
        names = [d or g for d, g in zip(ds.descriptions, GROUPS)]
        rw = from_bounds(*win.bounds, transform=ds.transform).round_offsets().round_lengths()
        arr = ds.read(window=rw, boundless=True, fill_value=0).astype("float32")
    if arr.shape[1:] != win.shape:
        raise ValueError(f"mask window {arr.shape[1:]} != {win.shape}")
    return {n: np.nan_to_num(arr[i]) for i, n in enumerate(names)}
