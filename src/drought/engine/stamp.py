"""Prove that the uploaded mask asset IS the local mask, then stamp it.

The GeoTIFF is uploaded by hand in the Code Editor (Assets → NEW → GeoTIFF;
name `grassland_mask_sg`; pyramiding policy SAMPLE). Uploading from Python would
need a Cloud Storage bucket, and with it billing.

So the check runs the other way round. A few hundred cells are drawn from the
local mask, with a fixed seed and weighted towards cells that hold grassland.
Their four purity values are then read back from the asset at the same cell
centres. Only if every value agrees (float32, to 1e-6) does the asset get the
local mask's content hash as its property `drought_mask_sha256`, which every
export checks before it starts.

A wrong file, a stale upload, a shifted grid, swapped bands or a 20 m mask all
fail. The comparison is a sample, not every pixel, but a mismatch anywhere
systematic shows up in a few hundred draws.
"""

from __future__ import annotations

import numpy as np

from ..config import CFG
from . import session
from .imagery import MASK_SHA_PROPERTY

__all__ = ["compare", "sample_cells", "stamp_mask"]

TOLERANCE = 1e-6


def sample_cells(path=None, n: int = 400, seed: int = 2026):
    """Cell centres (EPSG:2056) and their purity values, drawn from the local mask.

    Three quarters of the draws come from cells with any purity > 0, a quarter
    from anywhere, so both grassland and the zeros around it are compared.

    Returns:
        (xy array (n, 2), values array (n, bands), band names)
    """
    import rasterio

    from ..config import mask_path

    with rasterio.open(path or mask_path()) as ds:
        arr = ds.read()                               # bands, rows, cols
        tr = ds.transform
    rng = np.random.default_rng(seed)
    any_pos = np.flatnonzero(arr.max(axis=0).ravel() > 0)
    n_pos = min(len(any_pos), int(n * 0.75))
    idx = np.concatenate([rng.choice(any_pos, n_pos, replace=False),
                          rng.choice(arr[0].size, n - n_pos, replace=False)])
    rows, cols = np.unravel_index(idx, arr[0].shape)
    xs, ys = rasterio.transform.xy(tr, rows, cols, offset="center")
    vals = arr[:, rows, cols].T
    from ..mask.settings import GROUPS
    return np.column_stack([xs, ys]), vals, list(GROUPS)


def compare(local: np.ndarray, remote: dict[int, list[float]]) -> tuple[int, float]:
    """(cells missing from the asset, largest absolute difference)."""
    missing = [i for i in range(len(local)) if i not in remote]
    diffs = [np.max(np.abs(np.asarray(remote[i], float) - local[i]))
             for i in range(len(local)) if i in remote]
    return len(missing), float(max(diffs) if diffs else np.inf)


def stamp_mask(n: int = 400) -> str:
    """Compare the asset with the local mask; stamp it if they agree.

    Returns:
        The stamped hash.

    Raises:
        ProvenanceError: If the asset is missing, has the wrong bands, or differs.
    """
    from ..mask import verify_mask
    from ..provenance import ProvenanceError, mask_fingerprint, settings, settings_hash
    from .imagery import mask_image

    verify_mask()                                    # the local mask must be current
    fp = mask_fingerprint()
    name = CFG.gee.assets.mask
    info = session.asset_info(name)
    if info is None:
        raise ProvenanceError(name, "asset", None, "uploaded",
                              "Upload data/masks/grassland_mask_sg.tif in the Code Editor "
                              f"as an image asset named {name}, pyramiding policy SAMPLE.")
    bands = info.get("bands", [])
    if len(bands) != len(fp["bands"]):
        raise ProvenanceError(name, "bands", len(bands), len(fp["bands"]),
                              "The asset is not the four-band purity mask.")
    pyr = {b.get("pyramidingPolicy") for b in bands} - {None}
    if pyr and pyr != {"SAMPLE"}:
        raise ProvenanceError(name, "pyramidingPolicy", sorted(pyr), ["SAMPLE"],
                              "Re-upload with pyramiding policy SAMPLE: MEAN averages "
                              "purities at coarser scales.")

    xy, vals, band_names = sample_cells(n=n)
    ee = session.ee()
    proj = ee.Projection("EPSG:2056")
    fc = ee.FeatureCollection([ee.Feature(ee.Geometry.Point([float(x), float(y)], proj),
                                          {"i": i}) for i, (x, y) in enumerate(xy)])
    got = mask_image().sampleRegions(collection=fc, properties=["i"], scale=CFG.mask.res,
                                     projection=proj, geometries=False).getInfo()
    remote = {int(f["properties"]["i"]): [f["properties"].get(b, np.nan) for b in band_names]
              for f in got["features"]}
    missing, worst = compare(vals, remote)
    print(f"  compared {len(vals)} cells: {missing} missing, largest difference {worst:.2e}")
    if missing or not worst <= TOLERANCE:
        raise ProvenanceError(name, "pixel values", f"max |diff| {worst:.2e}, {missing} missing",
                              f"<= {TOLERANCE:g}", "The asset is not the local mask. "
                              "Upload the current GeoTIFF and run stamp-mask again.")
    ee.data.setAssetProperties(session.asset_id(name), {
        MASK_SHA_PROPERTY: fp["sha256"],
        "drought_mask_settings_hash": settings_hash(settings("mask")),
        "drought_mask_cells_compared": len(vals)})
    print(f"  stamped {session.asset_id(name)}: {MASK_SHA_PROPERTY} = {fp['sha256']}")
    return fp["sha256"]
