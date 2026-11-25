"""Per-parcel series for cut detection.

Mown parcels of at least cuts.min_area_ha (before the buffer), shrunk by
cuts.edge_buffer_m, measured on their pure mown cells: per stored mosaic in
cuts.doy, the median and count of NDVI and NDMI. Two sets: a random sample of
cuts.sample_n parcels (seeded) and every hand-labelled parcel. A mosaic whose
share of unusable grassland exceeds cuts.max_site_cloud is left out.

    <paths.series>/<site>/cuts_<year>.csv              the sample
    <paths.series>/<site>/cuts_<year>_labelled.csv     the parcels labelled that year
        parcel_id, area_ha, lnf_code, then <YYYYMMDDTHHMMSS>_SWISSEO_<index>_<stat>
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from drought.core.config import CFG, series_dir, years
from drought.products.common import cells, parcels, rasterize, write
from drought.sources import store

__all__ = ["labelled_ids", "output", "run"]

log = logging.getLogger(__name__)


def output(site: str, year: int, labelled: bool = False) -> Path:
    return series_dir(site) / f"cuts_{year}{'_labelled' if labelled else ''}.csv"


def labelled_ids(site: str) -> dict[int, list[str]]:
    """Hand-labelled parcel ids per year (empty without a label file)."""
    from drought.analysis.labels import load_labels
    out: dict[int, list[str]] = {}
    try:
        keys = load_labels(site)
    except FileNotFoundError:
        return out
    for key in keys:
        pid, year = key.rsplit("|", 1)
        out.setdefault(int(year), []).append(pid)
    return {y: sorted(set(v)) for y, v in out.items()}


def run(site: str, only_years: list[int] | None = None) -> list[Path]:
    """Compute the cut series of one site, per year."""
    c = cells(site)
    p = parcels(site, c.win)
    k = CFG.cuts
    mown = p[p.lnf_code.isin(list(CFG.landuse.mown)) & (p.area_ha >= float(k.min_area_ha))]
    mown = mown.sort_values("parcel_id").reset_index(drop=True)
    rng = np.random.default_rng(int(k.seed))
    n = min(int(k.sample_n), len(mown))
    sample = sorted(mown.parcel_id.iloc[rng.choice(len(mown), n, replace=False)]) if n else []
    labelled = labelled_ids(site)
    wanted = set(sample) | {i for v in labelled.values() for i in v}
    chosen = p[p.parcel_id.isin(wanted)].sort_values("parcel_id").reset_index(drop=True)
    missing = wanted - set(chosen.parcel_id)
    if missing:
        log.warning("  [!] %s: %d labelled parcel id(s) not in the cadastre inside the "
                    "window, e.g. %s", site, len(missing), sorted(missing)[:3])
    lab = rasterize(chosen.geometry.buffer(float(k.edge_buffer_m)),
                    np.arange(1, len(chosen) + 1), c.win)
    lab[c.purity["mown"] < float(CFG.mask.purity_min)] = 0
    flat = lab.ravel()
    has = flat > 0
    keys = flat[has]
    table = pd.DataFrame({"parcel_id": chosen.parcel_id, "area_ha": chosen.area_ha,
                          "lnf_code": chosen.lnf_code.astype("Int64")})
    log.info("  %s: %d sampled + labelled parcels -> %d", site, len(sample), len(chosen))

    out = []
    for year in only_years or years():
        cols: dict[str, np.ndarray] = {}
        for e, ix in store.iterate(site, year, c.win, doy=(int(k.doy[0]), int(k.doy[1]))):
            if c.cloud_pct(ix) > float(k.max_site_cloud):
                continue
            d = pd.DataFrame({"k": keys, "NDVI": ix["NDVI"].ravel()[has],
                              "NDMI": ix["NDMI"].ravel()[has]})
            g = d.groupby("k").agg(NDVI_median=("NDVI", "median"), NDVI_count=("NDVI", "count"),
                                   NDMI_median=("NDMI", "median"), NDMI_count=("NDMI", "count"))
            g = g.reindex(np.arange(1, len(chosen) + 1))
            tag = f"{e.datetime:%Y%m%dT%H%M%S}_SWISSEO"
            for col in g.columns:
                cols[f"{tag}_{col}"] = g[col].to_numpy()
        wide = pd.concat([table, pd.DataFrame(cols)], axis=1)
        out.append(write(wide[wide.parcel_id.isin(sample)], output(site, year), "cuts"))
        if labelled.get(year):
            out.append(write(wide[wide.parcel_id.isin(labelled[year])],
                             output(site, year, labelled=True), "cuts"))
    return out
