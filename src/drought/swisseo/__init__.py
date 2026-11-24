"""swissEO, read locally: the default source of the Sentinel-2 series.

swisstopo's swissEO S2-SR v200 is Sentinel-2 L2A surface reflectance for
Switzerland: one mosaic per orbit and date, co-registered to SWISSIMAGE,
reprojected to LV95 on a 10 m grid that coincides with the mask grid, with
OmniCloudMask cloud and shadow classes, a terrain-shadow layer and the scene
classification. Every band is a Cloud Optimised GeoTIFF, so only the site
windows are read, over HTTP. No account, no quota.

    drought swisseo items        what the catalogue holds for a site and year
    drought swisseo run          trajectory, sensitivity and cut series, one pass
    drought swisseo vhi          swisstopo's Vegetation Health Index per parcel

The outputs have the same columns and file names as the Earth Engine exports,
so everything downstream is unchanged. They land in paths.swisseo and are
registered with their provenance as they are written. Per-date results are
cached, so an interrupted run resumes where it stopped.
"""

from __future__ import annotations

__all__ = ["extract", "grid", "scene", "stac", "vhi"]
