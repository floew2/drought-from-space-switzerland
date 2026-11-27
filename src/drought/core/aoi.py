"""Areas of interest: the sites in config.yaml, or any polygon in a GeoJSON file.

    aoi.get("linthebene")                 a site from config.yaml
    aoi.resolve("my_area.geojson")        a polygon from a file: its name is
                                          the file name ("my_area")

A GeoJSON AOI is added to the configuration's sites for the rest of the run,
so every stage treats it exactly like a configured site.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from drought.core.config import CFG, Config

__all__ = ["AOI", "get", "resolve"]

_GEOMETRY: dict[str, Any] = {}       # name -> shapely geometry in LV95, for file AOIs


@dataclass(frozen=True)
class AOI:
    name: str
    label: str
    geometry: Any                       # shapely (Multi)Polygon in LV95 (EPSG:2056)
    stratum: str = "valley"

    @property
    def bounds(self) -> tuple[float, float, float, float]:
        return tuple(self.geometry.bounds)


def get(name: str) -> AOI:
    """A configured or registered AOI by name."""
    from pyproj import Transformer
    from shapely.geometry import box

    from drought.core.config import site_cfg

    s = site_cfg(name)
    geom = _GEOMETRY.get(name)
    if geom is None:
        tr = Transformer.from_crs(4326, 2056, always_xy=True)
        geom = box(*tr.transform_bounds(*s.bbox_wgs84))
    return AOI(name, s.label, geom, getattr(s, "stratum", "valley"))


def resolve(name_or_path: str) -> str:
    """A site name, or a GeoJSON file to register as one. Returns the name.

    Raises:
        KeyError: If it is neither a configured site nor an existing file.
        ValueError: If the file holds no polygon, or lies outside Switzerland.
    """
    p = Path(name_or_path)
    if p.suffix.lower() not in (".geojson", ".json") or not p.exists():
        get(name_or_path)                       # raises KeyError for an unknown name
        return name_or_path
    import geopandas as gpd

    g = gpd.read_file(p)
    if g.empty or not g.geom_type.isin(["Polygon", "MultiPolygon"]).any():
        raise ValueError(f"{p}: no polygon")
    if g.crs is None:
        g = g.set_crs(4326)                      # GeoJSON's default
    geom = g.to_crs(2056).union_all()
    w, s, e, n = g.to_crs(4326).total_bounds
    if not (5.9 <= w < e <= 10.6 and 45.8 <= s < n <= 47.9):
        raise ValueError(f"{p}: outside Switzerland; swissEO covers Switzerland only")
    name = p.stem
    setattr(CFG.sites, name, Config({"label": name, "note": f"from {p.name}",
                                     "bbox_wgs84": [float(w), float(s), float(e), float(n)],
                                     "stratum": "valley"}))
    _GEOMETRY[name] = geom
    return name
