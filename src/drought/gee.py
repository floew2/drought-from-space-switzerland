"""The CONFIG_HASH: one hash of every Earth Engine setting plus the local mask.

The Earth Engine exports (trajectory, purity sensitivity, cut series) run from
Python, in `drought.engine`. Each puts this hash in its file names
(`..._cfg<hash>`), so registering a download can tell whether it ran under the
current settings and mask. Every hash produced is kept in
`output/gee/manifest.json` with its settings, so a mismatch can name what
changed.

The hashed values are kept exactly as they were when interactive Code Editor
scripts still lived in gee/, so exports registered then keep their hash.
"""

from __future__ import annotations

import json
from pathlib import Path

from .config import CFG, ROOT

__all__ = ["expected_hash", "manifest_settings", "render_settings"]

# Settings of the Python exports that enter the CONFIG_HASH.
ENGINE_KEYS = ("gee.clear_threshold", "gee.strict_purity", "gee.parcels_year",
               "gee.trajectory", "gee.sensitivity", "gee.cuts")



def _shared() -> dict:
    """Method settings shared by every export."""
    from .mask.settings import GROUPS

    lu = CFG.landuse
    base = list(CFG.project.baseline_years)
    cur = int(CFG.project.current_year)
    return {
        "PURITY_MIN": CFG.mask.purity_min,
        "STRICT_PURITY": CFG.gee.strict_purity,
        "MASK_BANDS": list(GROUPS),
        "MOWN_CODES": list(lu.mown),
        "GRAZED_CODES": list(lu.grazed),
        "CODES": {"intensive": list(lu.intensive), "extensive": list(lu.extensive),
                  "mown": list(lu.mown), "grazed": list(lu.grazed)},
        "REGIMES": list(CFG.analysis.regimes),
        "YEARS": base + [cur],
        "YEAR_DRY": cur, "YEAR_CUR": cur,
        "BASELINE_YEARS": base,
    }


def render_settings(with_mask: bool = True) -> dict:
    """Every hashed setting, plus the mask fingerprint.

    This is what CONFIG_HASH hashes. Paths and the project name are excluded:
    they say where the assets live, not what the method is.
    """
    from .provenance import mask_fingerprint

    out = {k: v for k, v in _shared().items()}
    for key in ENGINE_KEYS:
        v = CFG.get(key)
        out[key] = v.to_dict() if hasattr(v, "to_dict") else v
    fp = mask_fingerprint() if with_mask else None
    out["mask.sha256"] = fp["sha256"] if fp else None
    return out


def expected_hash(with_mask: bool = True) -> str:
    """CONFIG_HASH for the current settings and mask."""
    from .provenance import settings_hash
    return settings_hash(render_settings(with_mask))


def _manifest_path(dest: Path) -> Path:
    return Path(dest) / "manifest.json"


def _remember(dest: Path, h: str) -> None:
    p = _manifest_path(dest)
    m = json.loads(p.read_text()) if p.exists() else {}
    m[h] = render_settings(with_mask=True)
    p.write_text(json.dumps(m, indent=1, sort_keys=True))


def manifest_settings(h: str, dest: str | Path | None = None) -> dict | None:
    """The settings a CONFIG_HASH was rendered with, if this machine rendered it."""
    for d in ([Path(dest)] if dest else [ROOT / "output" / "gee", ROOT / "output" / "synthetic" / "gee"]):
        p = _manifest_path(d)
        if p.exists():
            m = json.loads(p.read_text())
            if h in m:
                return m[h]
    return None
