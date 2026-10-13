"""Per-parcel NDVI / NDMI series: the input to cut detection.

Port of gee_cuts_extract.js, value for value. Per site, a seeded random sample
of mown parcels (at least gee.cuts.min_area_ha, shrunk by the edge buffer),
with their elevation; per year one export of the median and pixel count of
NDVI and NDMI over each parcel's pure cells, for every clear acquisition.

Why one image per year: the cost of the earlier designs was per-CALL overhead,
not pixel work (one reduceRegions per granule: ~8,500 EECU-s for 60 parcels and
one year; mosaicking by date first: ~31,800, because the computed mosaic is
not cached). So the season is stacked with toBands() into ONE image, and each
year costs ONE reduceRegions. The output is wide: <date>_<tile>_<index>_<stat>,
reshaped by `io.reshape_wide`.

Unlike the Code Editor script, the export names its columns explicitly. The
band names are read from Earth Engine first, so the CSV carries the parcel
attributes the pipeline uses and the series, and no geometry column (a cut
export with geometry must not enter the repository; provenance refuses it).
"""

from __future__ import annotations

from ..config import CFG
from ..io import INDEX_COLS, gee_dir
from . import imagery, session
from .tasks import Export, config_hash

__all__ = ["build", "exports", "prefix", "selectors", "years"]

STATS = ("median", "count")


def prefix(site: str, year: int, h: str, smoke: bool = False) -> str:
    return f"cuts_{site}_{year}{'_SMOKE' if smoke else ''}_cfg{h}"


def years() -> list[int]:
    y = [*CFG.project.baseline_years, CFG.project.current_year]
    return sorted(y, reverse=bool(CFG.gee.cuts.years_newest_first))


def selectors(band_names: list[str]) -> list[str]:
    """Columns of the export: parcel attributes, then median and count per band."""
    return [*INDEX_COLS, *(f"{b}_{s}" for b in band_names for s in STATS)]


def sample(site: str, n: int | None = None):
    """The parcel sample for a site, with area and elevation."""
    ee = session.ee()
    c = CFG.gee.cuts
    aoi = imagery.site_geometry(site)
    parcels = (ee.FeatureCollection(
        session.asset_id(f"{CFG.gee.assets.parcels}{CFG.gee.parcels_year}"))
        .filterBounds(aoi)
        .filter(ee.Filter.inList(CFG.gee.parcel_code_field, list(CFG.landuse.mown)))
        .map(lambda f: f.set("area_ha", f.geometry().area(1).divide(10000)))
        .filter(ee.Filter.gte("area_ha", c.min_area_ha))
        .map(lambda f: f.setGeometry(f.geometry().buffer(c.edge_buffer_m, 1))))
    # Cutting regime depends strongly on altitude, so every parcel carries it.
    elev = ee.Image("CGIAR/SRTM90_V4").select("elevation")
    with_elev = (elev.reduceRegions(collection=parcels,
                                    reducer=ee.Reducer.mean().setOutputs(["elev_m"]),
                                    scale=90, tileScale=4)
                 .filter(ee.Filter.notNull(["elev_m"])))
    code = CFG.gee.parcel_code_field
    return (with_elev.randomColumn("rnd", c.seed).sort("rnd").limit(n or c.sample_n)
            .map(lambda f: f.set("parcel_id", f.get("system:index"))
                 .set("lnf_code", f.get(code))))


def stack(site: str, year: int):
    """The season as one image: bands <date>_<tile>_<NDVI|NDMI>."""
    ee = session.ee()
    c = CFG.gee.cuts
    clear = CFG.gee.clear_threshold.cuts
    pure = imagery.mask_image().select("mown").gte(CFG.mask.purity_min)
    col = imagery.s2_linked(imagery.site_geometry(site), year, tuple(c.doy),
                            max_scene_cloud=c.max_scene_cloud, mgrs_tile=c.mgrs_tile)

    def bands(img):
        img = ee.Image(img)
        tag = img.date().format("YYYYMMdd").cat("_").cat(ee.String(img.get("MGRS_TILE")))
        return (imagery.index_bands(img.updateMask(pure), clear, ("NDVI", "NDMI"))
                .set("system:index", tag))

    return col.map(bands).toBands()


def build(site: str, year: int, n: int | None = None):
    """(FeatureCollection, selectors) for one site and year."""
    ee = session.ee()
    img = stack(site, year)
    names = img.bandNames().getInfo()          # metadata only: cheap
    rows = img.reduceRegions(collection=sample(site, n),
                             reducer=ee.Reducer.median().combine(ee.Reducer.count(), "", True),
                             scale=CFG.gee.cuts.scale_m, tileScale=8)
    return rows, selectors(names)


def exports(sites=None, only_years=None, smoke: bool = False) -> list[Export]:
    """One export per site and year. A smoke run: the current year, 100 parcels."""
    from ..config import ROOT

    h = config_hash()
    yrs = [CFG.project.current_year] if smoke else (only_years or years())
    out = []
    for site in sites or list(CFG.gee.cuts.sites):
        for year in yrs:
            fc, sel = build(site, year, n=100 if smoke else None)
            dest = ROOT / "output" / "ee_smoke" if smoke else gee_dir(site)
            out.append(Export(fc, prefix(site, year, h, smoke), dest,
                              None if smoke else "cuts", sel))
    return out
