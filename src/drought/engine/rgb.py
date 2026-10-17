"""A false-colour Sentinel-2 composite per site, for `drought mask overview`.

The first column of the imagery-to-mask figure. A cloud-masked median of one
summer (gee.rgb), bands NIR / red / green, on the 10 m LV95 grid of the mask,
exported as a GeoTIFF to Drive and downloaded to data/gee/<site>/rgb_<site>.tif.
For looking only: nothing in the analysis reads it, so it carries no
provenance record and does not enter the CONFIG_HASH.
"""

from __future__ import annotations

import math
from pathlib import Path

from ..config import CFG
from ..io import gee_dir
from . import imagery, session

__all__ = ["build", "run"]


def _grid(site: str) -> tuple[list[float], list[float]]:
    """Region corners and the 10 m transform, snapped to the mask grid."""
    from ..mask.io import site_bounds_lv95

    x0, y0, x1, y1 = site_bounds_lv95(site)
    x0, y0 = math.floor(x0 / 10) * 10, math.floor(y0 / 10) * 10
    x1, y1 = math.ceil(x1 / 10) * 10, math.ceil(y1 / 10) * 10
    return [x0, y0, x1, y1], [10, 0, x0, 0, -10, y1]


def build(site: str):
    """The composite for one site: median of clear pixels, bands NIR, red, green."""
    ee = session.ee()
    c = CFG.gee.rgb
    clear = CFG.gee.clear_threshold.map
    geom = imagery.site_geometry(site)
    col = imagery.s2_linked(geom, int(c.year), tuple(c.doy), max_scene_cloud=80)
    col = col.map(lambda i: ee.Image(i).updateMask(ee.Image(i).select("cs_cdf").gte(clear)))
    return col.select(list(c.bands)).median().toUint16()


def run(sites=None) -> list[Path]:
    """Export, wait, download. Returns the local GeoTIFFs."""
    from .tasks import Export, wait

    ee = session.ee()
    out, exps = [], []
    for site in sites or list(CFG.gee.trajectory.sites):
        dest = gee_dir(site) / f"rgb_{site}.tif"
        if dest.exists():
            print(f"  {dest.name} already here; skipped")
            out.append(dest)
            continue
        (x0, y0, x1, y1), transform = _grid(site)
        region = ee.Geometry.Rectangle([x0, y0, x1, y1], "EPSG:2056", False)
        name = f"rgb_{site}"
        task = ee.batch.Export.image.toDrive(
            image=build(site), description=name, folder=str(CFG.gee.drive_folder),
            fileNamePrefix=name, region=region, crs="EPSG:2056",
            crsTransform=transform, maxPixels=1e9, fileFormat="GeoTIFF")
        task.start()
        print(f"  started {name}  (task {task.id})")
        exps.append((Export(None, name, dest.parent, None, task=task), dest))
    wait([e for e, _ in exps])
    for e, dest in exps:
        found = session.drive_find(f"{e.prefix}.tif")
        if not found:
            raise FileNotFoundError(f"{e.prefix}.tif not found in Google Drive")
        session.drive_download(found[0]["id"], dest)
        print(f"  -> {dest}")
        out.append(dest)
    return out
