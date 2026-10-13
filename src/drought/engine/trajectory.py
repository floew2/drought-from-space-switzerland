"""Site-mean index trajectories: one row per site x regime x year x date.

Port of gee_grassland_trajectory.js, value for value. The mean of four indices
over the cells of one management system whose purity reaches mask.purity_min,
per clear acquisition, full season. One small CSV for all sites.

Cost: the mask is a stored asset and simply read. An early Code Editor version
rasterised parcels inside the per-image loop and burned ~97,000 EECU-seconds in
four minutes; reading the asset, the full run was ~4,000 (smoke test: 338).
"""

from __future__ import annotations

from ..config import CFG, ROOT
from . import imagery, session
from .tasks import Export, config_hash

__all__ = ["SELECTORS", "build", "exports", "prefix"]

SELECTORS = ["site", "regime", "year", "date", "doy", "NDVI", "kNDVI", "NDMI", "NDRE", "n_px"]
INDICES = ("NDVI", "kNDVI", "NDMI", "NDRE")


def prefix(h: str, smoke: bool = False) -> str:
    """Export file name (without .csv)."""
    p = round(CFG.mask.purity_min * 100)
    return f"grassland_trajectory_multi_p{p}{'_SMOKE' if smoke else ''}_cfg{h}"


def _years(smoke: bool) -> list[int]:
    years = [*CFG.project.baseline_years, CFG.project.current_year]
    return years[-2:] if smoke else years


def build(sites=None, years=None):
    """The FeatureCollection: one feature per site, regime, year and clear date."""
    ee = session.ee()
    c = CFG.gee.trajectory
    clear = CFG.gee.clear_threshold.trajectory
    mask = imagery.mask_image()
    reducer = ee.Reducer.mean().combine(ee.Reducer.count(), "", True)

    per_site = []
    for site in sites or list(c.sites):
        geom = imagery.site_geometry(site)
        for regime in CFG.analysis.regimes:
            cells = mask.select(regime).gte(CFG.mask.purity_min).selfMask()
            for year in years or _years(False):
                col = imagery.s2_linked(geom, year, tuple(c.doy))

                def row(img, site=site, regime=regime, year=year, cells=cells, geom=geom):
                    img = ee.Image(img)
                    idx = imagery.index_bands(img, clear, INDICES).updateMask(cells)
                    # ONE reduceRegion: the mean of all four bands plus a count.
                    st = idx.reduceRegion(reducer=reducer, geometry=geom,
                                          scale=c.scale_m, maxPixels=1e8,
                                          bestEffort=True, tileScale=4)
                    return ee.Feature(None, {
                        "site": site, "regime": regime, "year": year,
                        "date": img.date().format("YYYY-MM-dd"),
                        "doy": img.date().getRelative("day", "year").add(1),
                        **{n: st.get(f"{n}_mean") for n in INDICES},
                        # all four share one cloud mask, so one count serves
                        "n_px": st.get("NDVI_count"),
                    })

                per_site.append(col.map(row))
    return (ee.FeatureCollection(per_site).flatten()
            .filter(ee.Filter.notNull(list(INDICES)))
            .filter(ee.Filter.gte("n_px", c.min_pixels)))


def exports(smoke: bool = False) -> list[Export]:
    """The one export. A smoke run covers one site and two years, unregistered."""
    h = config_hash()
    sites = [CFG.gee.trajectory.sites[0]] if smoke else None
    dest = ROOT / CFG.paths.gee if not smoke else ROOT / "output" / "ee_smoke"
    return [Export(build(sites, _years(smoke)), prefix(h, smoke), dest,
                   None if smoke else "trajectory", SELECTORS)]
