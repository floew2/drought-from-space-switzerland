"""A false-colour Sentinel-2 image per site, for the imagery-to-mask figure.

The clearest swissEO mosaic in swisseo.rgb (year, DOY window) over the site
window, as NIR, red and green (B8, B4, B3) in swissEO's integer reflectance.
For looking only: nothing is computed from it.

    <paths.scenes>/<site>/rgb_<year>.tif
"""

from __future__ import annotations

import datetime as dt
import logging
from pathlib import Path

import numpy as np

from drought.core import grid
from drought.core.config import CFG, scenes_dir
from drought.sources import scene, stac, store

__all__ = ["BANDS", "fetch", "path"]

log = logging.getLogger(__name__)
BANDS = ("b08_10m", "b04_10m", "b03_10m")


def path(site: str) -> Path:
    return scenes_dir(site) / f"rgb_{int(CFG.swisseo.rgb.year)}.tif"


def fetch(site: str, session=None) -> Path:
    """Find the clearest mosaic and store its false-colour image. Returns the file."""
    r = CFG.swisseo.rgb
    year = int(r.year)
    jan1 = dt.date(year, 1, 1)
    its = stac.search(getattr(CFG.sites, site).bbox_wgs84,
                      jan1 + dt.timedelta(days=int(r.doy[0]) - 1),
                      jan1 + dt.timedelta(days=int(r.doy[1]) - 1), session=session)
    if not its:
        raise FileNotFoundError(f"{site}: no mosaic in swisseo.rgb's window")
    win = grid.site_window(site)
    with scene.gdal_env():
        clear = [(float((scene.read_band(i.href(scene.CLOUD), win, fill=255) == 0).mean()), i)
                 for i in its]
    share, best = max(clear, key=lambda t: t[0])
    a = store._retry(lambda: np.stack([scene.read_band(best.href(b), win) for b in BANDS]))
    out = store._write(path(site), a.astype(np.uint16), win, ("B8", "B4", "B3"), 0,
                       {"datetime": best.datetime.isoformat(), "clear_share": f"{share:.3f}"})
    log.info("  %s: %s, %.0f%% clear -> %s", site, best.id, 100 * share, out)
    return out
