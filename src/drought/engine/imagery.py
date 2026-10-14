"""What the three exports share: the mask asset, sites, and Sentinel-2 indices.

The mask asset is checked before anything runs. `stamp.stamp_mask` compares the
uploaded asset with the local mask pixel by pixel and, if they agree, stores the
local mask's content hash as the asset property `drought_mask_sha256`. If the
asset carries another hash, or none, an export refuses to start. Without this,
a run could read an older mask while the local files, and the provenance
records, described a newer one.
"""

from __future__ import annotations

from ..config import CFG
from . import session

__all__ = ["MASK_SHA_PROPERTY", "check_mask_asset", "index_bands", "mask_image",
           "s2_linked", "site_geometry"]

MASK_SHA_PROPERTY = "drought_mask_sha256"
S2 = "COPERNICUS/S2_SR_HARMONIZED"
CS_PLUS = "GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED"


def check_mask_asset(properties: dict | None = None) -> str:
    """Refuse to run against a mask asset that is not the local mask.

    Args:
        properties: The asset's properties (default: read from Earth Engine).

    Returns:
        The mask hash both sides agree on.

    Raises:
        ProvenanceError: If the asset is missing, unstamped, or another build.
    """
    from ..provenance import ProvenanceError, mask_fingerprint

    fp = mask_fingerprint()
    name = CFG.gee.assets.mask
    if fp is None:
        raise ProvenanceError(name, "local mask", None, "a built mask",
                              "Build it first: drought mask build")
    if properties is None:
        info = session.asset_info(name)
        props = None if info is None else info.get("properties", {})
    else:
        props = properties
    if props is None:
        raise ProvenanceError(name, "asset", None, fp["sha256"],
                              "No mask asset in Earth Engine. Upload the GeoTIFF in the "
                              "Code Editor, then: drought ee stamp-mask")
    remote = props.get(MASK_SHA_PROPERTY)
    if remote != fp["sha256"]:
        raise ProvenanceError(name, MASK_SHA_PROPERTY, remote, fp["sha256"],
                              "The Earth Engine asset is not stamped as the local mask. "
                              "Upload the current GeoTIFF in the Code Editor, then: "
                              "drought ee stamp-mask")
    return remote


def mask_image():
    """The four purity bands, by name."""
    from ..mask.settings import GROUPS

    ee = session.ee()
    return ee.Image(session.asset_id(CFG.gee.assets.mask)).rename(list(GROUPS))


def site_geometry(site: str):
    """The site rectangle (WGS84) from config.yaml."""
    ee = session.ee()
    return ee.Geometry.Rectangle(list(getattr(CFG.sites, site).bbox_wgs84))


def s2_linked(geom, year: int, doy: tuple[int, int], max_scene_cloud: float | None = None,
              mgrs_tile: str | None = None):
    """Sentinel-2 L2A for one year and day-of-year window, with Cloud Score+ linked."""
    ee = session.ee()
    col = (ee.ImageCollection(S2).filterBounds(geom)
           .filter(ee.Filter.calendarRange(year, year, "year"))
           .filter(ee.Filter.calendarRange(int(doy[0]), int(doy[1]), "day_of_year")))
    if max_scene_cloud is not None:
        col = col.filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", max_scene_cloud))
    if mgrs_tile:
        col = col.filter(ee.Filter.eq("MGRS_TILE", str(mgrs_tile).lstrip("T")))
    return col.linkCollection(ee.ImageCollection(CS_PLUS), ["cs_cdf"])


def index_bands(img, clear: float, names=("NDVI", "kNDVI", "NDMI", "NDRE")):
    """Cloud-masked index bands of one linked image.

    kNDVI with sigma = 0.5 (NIR + Red) reduces to tanh(NDVI^2).
    """
    m = img.updateMask(img.select("cs_cdf").gte(clear))
    ndvi = m.normalizedDifference(["B8", "B4"]).rename("NDVI")
    bands = {"NDVI": ndvi,
             "kNDVI": ndvi.pow(2).tanh().rename("kNDVI"),
             "NDMI": m.normalizedDifference(["B8A", "B11"]).rename("NDMI"),
             "NDRE": m.normalizedDifference(["B8A", "B5"]).rename("NDRE")}
    out = bands[names[0]]
    for n in names[1:]:
        out = out.addBands(bands[n])
    return out
