"""Earth Engine, run from Python: the three exports, and a check of the mask asset.

These replace the Code Editor scripts gee_grassland_trajectory.js,
gee_purity_sensitivity.js and gee_cuts_extract.js. Every setting is read from
config.yaml, so nothing has to be pasted, rendered or kept in step by hand.

    drought ee check          project, Drive, and mask asset == local mask?
    drought ee stamp-mask     compare the uploaded asset with the local mask; stamp it
    drought ee trajectory     site-mean index series      -> data/gee/
    drought ee sensitivity    purity-class series          -> data/gee/
    drought ee cuts           per-parcel series, per year  -> data/gee/<site>/

Exports run as batch tasks to Google Drive and are downloaded from there with
the Earth Engine login, which already includes Drive access. No Cloud Storage
bucket, no billing.

Each export first checks that the mask asset in Earth Engine IS the local mask:
`stamp-mask` compares a few hundred cells of the asset with the local GeoTIFF
and, only if they agree, stores the mask's content hash on the asset; an
export refuses to start if the stamp differs. Each download is registered with
its provenance record at once, so a stale or foreign export can never reach
the analysis (METHOD §6, bug 6).

Needs the optional dependency: `uv sync --extra ee`, then once
`earthengine authenticate`. The project goes into config.local.yaml (not
committed): `gee: {project: <your-ee-project>}`.
"""

from __future__ import annotations

__all__ = ["cuts", "imagery", "sensitivity", "session", "stamp", "tasks", "trajectory"]
