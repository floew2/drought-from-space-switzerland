"""Purity sensitivity: does the anomaly depend on where the threshold is set?

Port of gee_purity_sensitivity.js, value for value.

One extraction, not one per threshold. Purity is quantised into classes, and a
grouped reducer returns the mean index and pixel count per class for each
year's median composite. The mean at any threshold t is then the
pixel-weighted mean of the classes at or above t, reconstructed by
`figures.figure_purity_sensitivity`, identical to re-extracting at t.

Which parcels enter: those holding at least one cell at gee.strict_purity. The
default 0.0 admits every parcel, so the sweep covers the population the anomaly
is reported for. If the curve drifts, a run at 0.90 separates contamination
from composition (METHOD §2.5).
"""

from __future__ import annotations

from ..config import CFG, ROOT
from . import imagery, session
from .tasks import Export, config_hash

__all__ = ["SELECTORS", "build", "exports", "prefix"]

SELECTORS = ["year", "pclass", "purity_lo", "mean", "n_px"]


def prefix(site: str, regime: str, h: str) -> str:
    return f"purity_sensitivity_{site}_{regime}_cfg{h}"


def build(site: str, regime: str):
    """The FeatureCollection: one feature per year and purity class."""
    ee = session.ee()
    c = CFG.gee.sensitivity
    clear = CFG.gee.clear_threshold.sensitivity
    lu = CFG.landuse
    doy = tuple(getattr(c.windows, site))
    aoi = imagery.site_geometry(site)
    purity = imagery.mask_image().select(regime)

    parcels = (ee.FeatureCollection(
        session.asset_id(f"{CFG.gee.assets.parcels}{CFG.gee.parcels_year}"))
        .filterBounds(aoi)
        .filter(ee.Filter.inList(CFG.gee.parcel_code_field, list(getattr(lu, regime)))))

    # A parcel qualifies if ANY cell inside it reaches the strict purity.
    qualified = (purity.gte(CFG.gee.strict_purity)
                 .reduceRegions(collection=parcels,
                                reducer=ee.Reducer.max().setOutputs(["pure"]),
                                scale=c.scale_m, tileScale=4)
                 .filter(ee.Filter.eq("pure", 1)))
    eligible = ee.Image(0).byte().paint(qualified, 1).selfMask()

    n = int(c.n_classes)
    pclass = (purity.multiply(n).floor().min(n - 1)
              .updateMask(purity.gt(0)).updateMask(eligible)
              .rename("pclass").toInt())

    def per_year(year):
        col = imagery.s2_linked(aoi, year, doy, max_scene_cloud=c.max_scene_cloud)
        col = col.map(lambda img: imagery.index_bands(ee.Image(img), clear, (c.index,))
                      .rename("v"))
        comp = col.median().rename("v").updateMask(pclass.mask())
        stats = comp.addBands(pclass).reduceRegion(
            reducer=ee.Reducer.mean().combine(ee.Reducer.count(), "", True)
            .group(groupField=1, groupName="pclass"),
            geometry=aoi, scale=c.scale_m, maxPixels=1e9, bestEffort=True, tileScale=4)

        def feat(g):
            g = ee.Dictionary(g)
            return ee.Feature(None, {"year": year, "pclass": g.get("pclass"),
                                     "purity_lo": ee.Number(g.get("pclass")).divide(n),
                                     "mean": g.get("mean"), "n_px": g.get("count")})

        return ee.FeatureCollection(ee.List(stats.get("groups")).map(feat))

    years = [CFG.project.current_year, *CFG.project.baseline_years]
    return ee.FeatureCollection([per_year(y) for y in years]).flatten()


def exports(sites=None, regimes=None) -> list[Export]:
    """One export per site and management system."""
    h = config_hash()
    out = []
    for site in sites or list(CFG.gee.sensitivity.windows.to_dict()):
        for regime in regimes or list(CFG.analysis.regimes):
            out.append(Export(build(site, regime), prefix(site, regime, h),
                              ROOT / CFG.paths.gee, "sensitivity", SELECTORS))
    return out
