"""A stable id for every parcel, stored in the parcel data.

The id is the parcel centroid in LV95, rounded to the metre: "2718345_1226789".
It depends on the geometry alone, so any tool on any platform gives the same id
for the same field, and pasting it into a GIS finds the field. Labels and every
per-parcel output refer to this id and nothing else.

    drought parcels ids     write data/parcels/parcels_<year>.shp with a parcel_id column

The written file is the one the per-parcel work reads, and the one to upload
to Earth Engine (asset gee.assets.parcels + year) when that backend is used.
"""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd

from .config import CFG, ROOT

__all__ = ["make_ids", "path", "read", "write"]


def make_ids(geoms: gpd.GeoSeries) -> list[str]:
    """Centroid in LV95, rounded to the metre, as "<E>_<N>"."""
    c = geoms.to_crs(2056).centroid
    return [f"{x:.0f}_{y:.0f}" for x, y in zip(c.x, c.y)]


def path(year: int | None = None) -> Path:
    year = int(year or CFG.gee.parcels_year)
    return ROOT / CFG.paths.parcel_ids / f"parcels_{year}.shp"


def write(year: int | None = None) -> Path:
    """Read the delivery for `year`, add parcel_id, write it. Fails on duplicates."""
    import warnings

    from .mask.parcels import _year_files
    from .mask.settings import CODE_FIELD_CANDIDATES

    year = int(year or CFG.gee.parcels_year)
    src = _year_files()[year]
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*organizePolygons|.*Simple Geometry")
        g = gpd.read_file(src)
    if g.crs is None:
        raise ValueError(f"{src}: no CRS")
    code = next((c for c in CODE_FIELD_CANDIDATES if c in g.columns), None)
    if code is None:
        raise ValueError(f"{src}: none of {CODE_FIELD_CANDIDATES} in {list(g.columns)}")
    g = g[g.geometry.notna() & ~g.geometry.is_empty].to_crs(2056)
    out = gpd.GeoDataFrame({"parcel_id": make_ids(g.geometry), "lnf_code": g[code].values},
                           geometry=g.geometry.values, crs=2056)
    dup = out.parcel_id[out.parcel_id.duplicated(keep=False)]
    if len(dup):
        raise ValueError(f"{len(dup)} parcels share a centroid id, e.g. "
                         f"{sorted(set(dup))[:5]}: check the delivery for duplicate geometries")
    p = path(year)
    p.parent.mkdir(parents=True, exist_ok=True)
    out.to_file(p)
    print(f"  {len(out):,} parcels from {src.name} -> {p}")
    return p


def read(year: int | None = None) -> gpd.GeoDataFrame:
    """The parcels with their ids. Run `drought parcels ids` first."""
    p = path(year)
    if not p.exists():
        raise SystemExit(f"{p} not found: run `drought parcels ids` first")
    g = gpd.read_file(p)
    g["parcel_id"] = g.parcel_id.astype(str)
    return g
