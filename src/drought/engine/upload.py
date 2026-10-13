"""Upload the local mask to Earth Engine, stamped with its content hash.

The GeoTIFF goes to the bucket, then into an image asset through an ingestion
manifest. That fixes three things the Code Editor upload left to hand:

- the band names (no renaming by position downstream),
- pyramidingPolicy SAMPLE for every band (the default for float, MEAN, would
  average purities at coarser scales, giving a purity of nothing),
- the asset property `drought_mask_sha256`, which every export checks against
  the local mask before it starts (`imagery.check_mask_asset`).
"""

from __future__ import annotations

import time

from ..config import CFG
from . import session
from .imagery import MASK_SHA_PROPERTY

__all__ = ["manifest", "upload_mask"]


def manifest(uri: str, fp: dict, settings_hash: str) -> dict:
    """The ingestion manifest for the mask asset."""
    from ..mask.settings import GROUPS

    return {
        "name": session.asset_id(CFG.gee.assets.mask),
        "tilesets": [{"id": "mask", "sources": [{"uris": [uri]}]}],
        "bands": [{"id": b, "tilesetId": "mask", "tilesetBandIndex": i,
                   "pyramidingPolicy": "SAMPLE"} for i, b in enumerate(GROUPS)],
        "properties": {MASK_SHA_PROPERTY: fp["sha256"],
                       "drought_mask_settings_hash": settings_hash,
                       "drought_mask_res_m": fp["res_m"]},
    }


def upload_mask(poll_s: int = 15) -> str:
    """Upload data/masks/grassland_mask_sg.tif and replace the asset.

    Returns:
        The asset id.
    """
    from ..config import mask_path
    from ..mask import verify_mask
    from ..provenance import mask_fingerprint, settings, settings_hash

    verify_mask()                                  # refuses a stale or foreign mask
    fp = mask_fingerprint()
    ee = session.ee()
    name = session.object_name("mask", f"grassland_mask_sg_{fp['sha256']}.tif")
    if name in session.list_objects(name):
        print(f"  gs://{session.bucket()}/{name} already uploaded")
        uri = f"gs://{session.bucket()}/{name}"
    else:
        print(f"  uploading {mask_path().name} ({mask_path().stat().st_size / 1e6:.1f} MB)")
        uri = session.upload(mask_path(), name, "image/tiff")
    m = manifest(uri, fp, settings_hash(settings("mask")))
    task_id = ee.data.newTaskId()[0]
    ee.data.startIngestion(task_id, m, allow_overwrite=True)
    print(f"  ingesting into {m['name']} (task {task_id})")
    while True:
        st = ee.data.getTaskStatus(task_id)[0]
        if st["state"] == "COMPLETED":
            break
        if st["state"] in ("FAILED", "CANCELLED"):
            raise RuntimeError(f"ingestion {st['state']}: {st.get('error_message', '')}")
        time.sleep(poll_s)
    print(f"  done: {m['name']}, {MASK_SHA_PROPERTY} = {fp['sha256']}")
    return m["name"]
