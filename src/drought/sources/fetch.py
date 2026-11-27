"""Download swissEO for an AOI and years into the store. Needs no mask, no parcels.

    fetch.run(["linthebene"], [2025, 2026])

For each site and year it asks the catalogue which mosaics exist, records the
answer (store.write_catalogue), and downloads the ones not yet on disk. A
mosaic that still fails after store.ATTEMPTS is reported and left missing; the
next run fetches only what is missing. The VHI maps are fetched for the
current year, over the site's window.
"""

from __future__ import annotations

import datetime as dt
import logging
from concurrent.futures import ThreadPoolExecutor

from drought.core import grid
from drought.core.config import CFG, window
from drought.sources import stac, store

__all__ = ["items", "run"]

log = logging.getLogger(__name__)


def _dates(year: int, doy: tuple[int, int]) -> tuple[dt.date, dt.date]:
    jan1 = dt.date(year, 1, 1)
    return jan1 + dt.timedelta(days=int(doy[0]) - 1), jan1 + dt.timedelta(days=int(doy[1]) - 1)


def items(site: str, year: int, kind: str = "s2", session=None) -> list[stac.Item]:
    """What the catalogue holds for a site and year (S2: the season; VHI: the window)."""
    bbox = getattr(CFG.sites, site).bbox_wgs84
    if kind == "s2":
        return stac.search(bbox, *_dates(year, CFG.swisseo.season_doy), session=session)
    return stac.search(bbox, *_dates(year, window(site)),
                       collection=CFG.swisseo.vhi.collection, session=session)


def _download(site: str, year: int, kind: str, its: list, workers: int) -> list[str]:
    todo = [it for it in its if not store.path(site, it.id, kind).exists()]
    if not todo:
        log.info("  %s %s %d: all %d on disk", site, kind, year, len(its))
        return []
    win = grid.site_window(site)
    get = store.fetch if kind == "s2" else store.fetch_vhi

    def one(it):
        try:
            get(it, win, site)
            return None
        except Exception as e:                  # noqa: BLE001 - reported, run continues
            log.warning("  [!] %s %s: %s: %s", site, it.id, type(e).__name__, e)
            return it.id

    with ThreadPoolExecutor(max_workers=max(1, workers)) as ex:
        bad = [r for r in ex.map(one, todo) if r]
    log.info("  %s %s %d: %d of %d downloaded (%d already on disk)", site, kind, year,
             len(todo) - len(bad), len(todo), len(its) - len(todo))
    return bad


def run(sites: list[str], years: list[int], workers: int | None = None,
        vhi: bool = True, dry_run: bool = False, session=None) -> list[str]:
    """Fetch S2 (and VHI) for every site and year.

    Args:
        sites: Site names (configured or registered, see core.aoi).
        years: Calendar years.
        workers: Mosaics downloaded in parallel (default swisseo.workers).
        vhi: Also fetch the VHI maps for the current year, if it is in `years`.
        dry_run: Count what the catalogue holds and what is on disk; download nothing.
        session: For tests.

    Returns:
        The ids that could not be downloaded (empty when everything is on disk).
    """
    workers = workers or int(CFG.swisseo.workers)
    cur = int(CFG.project.current_year)
    failed: list[str] = []
    for site in sites:
        jobs = [("s2", y) for y in years] + ([("vhi", cur)] if vhi and cur in years else [])
        for kind, year in jobs:
            its = items(site, year, kind, session)
            if dry_run:
                have = sum(store.path(site, i.id, kind).exists() for i in its)
                log.info("  %s %s %d: %d mosaics, %d on disk", site, kind, year, len(its), have)
                continue
            store.write_catalogue(site, year, its, kind)
            failed += _download(site, year, kind, its, workers)
    if failed:
        log.warning("  %d mosaic(s) could not be downloaded, e.g. %s. Run the same "
                    "command again to fetch them.", len(failed), failed[:3])
    return failed
