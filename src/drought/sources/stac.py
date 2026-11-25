"""The swisstopo STAC catalogue: which mosaics exist for a site and period."""

from __future__ import annotations

import datetime as dt
from dataclasses import dataclass, field

from drought.core.config import CFG

__all__ = ["Item", "search"]

TIMEOUT_S = 60


@dataclass(frozen=True)
class Item:
    """One mosaic: an orbit on a date. `assets` maps a short key to its URL."""
    id: str
    datetime: dt.datetime
    assets: dict = field(hash=False, compare=False)

    @property
    def doy(self) -> int:
        return self.datetime.timetuple().tm_yday

    def href(self, key: str) -> str:
        """URL of the asset whose name ends in `_<key>.tif` (e.g. "b08_10m")."""
        try:
            return self.assets[key]
        except KeyError:
            raise KeyError(f"{self.id}: no asset {key!r}; has {sorted(self.assets)}") from None


def _short(name: str) -> str:
    """swisseo_s2-sr_v200_mosaic_2026-07-02t103021_b08_10m.tif -> b08_10m."""
    stem = name.rsplit("/", 1)[-1].rsplit(".", 1)[0]
    parts = stem.split("_")
    return "_".join(parts[-2:]) if parts[-1].endswith("m") and parts[-1][:-1].isdigit() \
        else parts[-1]


def _parse(feature: dict) -> Item:
    t = dt.datetime.fromisoformat(feature["properties"]["datetime"].replace("Z", "+00:00"))
    assets = {_short(k): v["href"] for k, v in feature["assets"].items()}
    return Item(feature["id"], t, assets)


def search(bbox_wgs84, start: dt.date, end: dt.date, collection: str | None = None,
           session=None) -> list[Item]:
    """Every mosaic intersecting `bbox_wgs84` between `start` and `end` (inclusive).

    Args:
        bbox_wgs84: (west, south, east, north).
        start, end: Dates, inclusive.
        collection: Default: swisseo.collection.
        session: A requests-like session (tests pass a fake).

    Returns:
        Items sorted by time.
    """
    import requests

    s = session or requests.Session()
    coll = collection or CFG.swisseo.collection
    url: str | None = f"{CFG.swisseo.stac.rstrip('/')}/collections/{coll}/items"
    params: dict[str, str | int] | None = {"bbox": ",".join(f"{v:g}" for v in bbox_wgs84),
              "datetime": f"{start:%Y-%m-%d}T00:00:00Z/{end:%Y-%m-%d}T23:59:59Z",
              "limit": 100}
    out = []
    while url:
        r = s.get(url, params=params, timeout=TIMEOUT_S)
        r.raise_for_status()
        page = r.json()
        out += [_parse(f) for f in page.get("features", [])]
        nxt = [lk for lk in page.get("links", []) if lk.get("rel") == "next"]
        url, params = (nxt[0]["href"], None) if nxt else (None, None)
    return sorted(out, key=lambda i: (i.datetime, i.id))
