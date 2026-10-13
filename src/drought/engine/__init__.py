"""Earth Engine, run from Python: the three exports and the mask upload.

These replace the Code Editor scripts gee_grassland_trajectory.js,
gee_purity_sensitivity.js and gee_cuts_extract.js. Every setting is read from
config.yaml, so nothing has to be pasted, rendered or kept in step by hand.

    drought ee check          project, bucket and mask asset in order?
    drought ee upload-mask    local mask -> bucket -> Earth Engine asset
    drought ee trajectory     site-mean index series      -> data/gee/
    drought ee sensitivity    purity-class series          -> data/gee/
    drought ee cuts           per-parcel series, per year  -> data/gee/<site>/
    drought ee fetch          download finished exports of the current settings

Each export first checks that the mask asset in Earth Engine IS the local mask:
the upload stores the mask's content hash as an asset property, and an export
refuses to start if it differs. Each download is registered with its
provenance record at once, so a stale or foreign export can never reach the
analysis (METHOD §6, bug 6).

Needs the optional dependency: `uv sync --extra ee`, then once
`earthengine authenticate`. The project and bucket go into config.local.yaml
(not committed):

    gee:
      project: <your-ee-project>
      bucket: <a-bucket-in-that-project>
"""

from __future__ import annotations

__all__ = ["cuts", "imagery", "sensitivity", "session", "tasks", "trajectory", "upload"]
