"""Render Earth Engine settings from config.yaml into the scripts in gee/.

Every setting the scripts share with the Python side — the purity threshold,
the mask band order, the land-use code lists, the years, the Cloud Score+
threshold — is written into a marked block at the top of each script. The
scripts no longer define these values themselves, so they cannot drift from
config.yaml. CI runs `drought gee render --check` and fails if a committed
script differs from what the renderer produces.

Two kinds of output:

gee/*.js            committed. The Earth Engine project is a placeholder and no
                    mask is assumed, so the files carry nothing private and are
                    reproducible anywhere.
<out>/*.js          for pasting into the Code Editor: `drought --set
                    gee.project=<project> gee render --out output/gee`. Carries
                    the real project and the fingerprint of the local mask.

Each rendered copy has a CONFIG_HASH: a hash of every rendered setting plus the
mask fingerprint. The data exports put it in their file names, so registering
a download (`drought provenance register`) can tell whether the export ran
under the current settings and mask. Every hash rendered is kept in
`<out>/manifest.json` with its settings, so a mismatch can name what changed.
"""

from __future__ import annotations

import json
from pathlib import Path

from .config import CFG, ROOT

__all__ = ["SCRIPTS", "check", "expected_hash", "manifest_settings", "render",
           "render_settings"]

GEE_DIR = ROOT / "gee"
BEGIN = "// >>> drought gee render"
END = "// <<< end of generated block"
PLACEHOLDER = "YOUR_EE_PROJECT"

# Which generated variables each script receives. gee_switzerland_rgb.js is not
# rendered: it is MODIS, with a deliberately different baseline (2011-2020).
SCRIPTS = {
    "gee_cuts_extract.js": ("PROJECT", "ASSETS", "CONFIG_HASH", "PURITY_MIN",
                            "CLEAR_THRESHOLD", "MASK_BANDS", "MOWN_CODES"),
    "gee_grassland_map.js": ("PROJECT", "ASSETS", "CONFIG_HASH", "PURITY_MIN",
                             "CLEAR_THRESHOLD", "MASK_BANDS", "MOWN_CODES", "GRAZED_CODES",
                             "YEAR_DRY", "BASELINE_YEARS"),
    "gee_grassland_trajectory.js": ("PROJECT", "ASSETS", "CONFIG_HASH", "PURITY_MIN",
                                    "CLEAR_THRESHOLD", "MASK_BANDS", "REGIMES", "YEARS"),
    "gee_purity_sensitivity.js": ("PROJECT", "ASSETS", "CONFIG_HASH", "PURITY_MIN",
                                  "CLEAR_THRESHOLD", "MASK_BANDS", "CODES", "YEAR_CUR",
                                  "BASELINE_YEARS"),
    "gee_s2_hotspots.js": ("CLEAR_THRESHOLD", "YEAR_DRY", "BASELINE_YEARS"),
}

COMMENTS = {
    "PROJECT": "Earth Engine project with the assets",
    "CONFIG_HASH": "settings + mask this copy was rendered with; goes into export names",
    "PURITY_MIN": "mask.purity_min",
    "CLEAR_THRESHOLD": "Cloud Score+ cs_cdf, gee.clear_threshold",
    "MASK_BANDS": "band order of the mask raster, as written by build_mask()",
    "MOWN_CODES": "landuse.mown", "GRAZED_CODES": "landuse.grazed",
    "CODES": "landuse", "REGIMES": "analysis.regimes",
    "YEARS": "project.baseline_years + current_year",
    "YEAR_DRY": "project.current_year", "YEAR_CUR": "project.current_year",
    "BASELINE_YEARS": "project.baseline_years",
}


def _shared() -> dict:
    """Generated values that do not depend on the script."""
    from .mask.settings import GROUPS

    lu = CFG.landuse
    base = list(CFG.project.baseline_years)
    cur = int(CFG.project.current_year)
    return {
        "PURITY_MIN": CFG.mask.purity_min,
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
    """Every rendered value across all scripts, plus the mask fingerprint.

    This is what CONFIG_HASH hashes. Paths and the project name are excluded:
    they say where the assets live, not what the method is.
    """
    from .provenance import mask_fingerprint

    out = {k: v for k, v in _shared().items()}
    for script in SCRIPTS:
        out[f"CLEAR_THRESHOLD.{Path(script).stem}"] = getattr(
            CFG.gee.clear_threshold, Path(script).stem)
    fp = mask_fingerprint() if with_mask else None
    out["mask.sha256"] = fp["sha256"] if fp else None
    return out


def expected_hash(with_mask: bool = True) -> str:
    """CONFIG_HASH for the current settings and mask."""
    from .provenance import settings_hash
    return settings_hash(render_settings(with_mask))


def _values(script: str, project: str | None, config_hash: str) -> dict:
    v = _shared()
    v["PROJECT"] = project or PLACEHOLDER
    v["CONFIG_HASH"] = config_hash
    v["CLEAR_THRESHOLD"] = getattr(CFG.gee.clear_threshold, Path(script).stem)
    return {k: v[k] for k in SCRIPTS[script] if k != "ASSETS"}


def block(script: str, project: str | None = None, config_hash: str = "") -> str:
    """The generated block for one script."""
    lines = [f"{BEGIN}: generated from config.yaml. Do not edit by hand —",
             "// change config.yaml and run `drought gee render`.",
             "// Paste-ready copy: drought --set gee.project=<project> gee render --out output/gee"]
    for name, value in _values(script, project, config_hash).items():
        js = json.dumps(value)
        if name in ("CODES",):
            js = "{\n" + ",\n".join(f"  {k}: {json.dumps(v)}" for k, v in value.items()) + "\n}"
        lines.append(f"var {name} = {js};" + (f"   // {COMMENTS[name]}" if name in COMMENTS
                                               and "\n" not in js else ""))
        if name == "PROJECT" and "ASSETS" in SCRIPTS[script]:
            lines.append("var ASSETS = 'projects/' + PROJECT + '/assets/';")
            lines.append(f"var MASK_ASSET = ASSETS + {json.dumps(CFG.gee.assets.mask)};"
                         "   // uploaded from build_mask()'s GeoTIFF")
            lines.append(f"var PARCEL_PREFIX = ASSETS + {json.dumps(CFG.gee.assets.parcels)};"
                         "   // + year")
    lines.append(END)
    return "\n".join(lines)


def apply(text: str, blk: str) -> str:
    """Put the block into a script: replace an existing one, or insert it
    after the header comment."""
    if BEGIN in text:
        a = text.index(BEGIN)
        b = text.index(END, a) + len(END)
        return text[:a] + blk + text[b:]
    lines = text.split("\n")
    i = next(k for k, ln in enumerate(lines) if ln.strip() and not ln.lstrip().startswith("//"))
    return "\n".join(lines[:i] + [blk, ""] + lines[i:])


def render(out: str | Path | None = None, check_only: bool = False) -> list[Path]:
    """Render every script.

    Args:
        out: Write paste-ready copies here, with the configured project and
            the local mask's fingerprint. None updates the committed templates
            in gee/ (placeholder project, no mask).
        check_only: Compare the templates in gee/ with a fresh render and
            write nothing.

    Returns:
        Paths written, or (with check_only) paths that differ.

    Raises:
        ValueError: If paste-ready copies are requested without a project.
    """
    template = out is None
    project = None if template else CFG.gee.project
    if not template and not project:
        raise ValueError("set the project: drought --set gee.project=<project> gee render "
                         "--out output/gee")
    h = expected_hash(with_mask=not template)
    dest = GEE_DIR if template else Path(out)
    dest.mkdir(parents=True, exist_ok=True)
    changed = []
    for script in SCRIPTS:
        src = (GEE_DIR / script).read_text()
        new = apply(src, block(script, project, h))
        target = dest / script
        if check_only:
            if new != src:
                changed.append(target)
            continue
        if not target.exists() or target.read_text() != new:
            target.write_text(new)
            changed.append(target)
    if not check_only and not template:
        _remember(dest, h)
    return changed


def check() -> list[Path]:
    """Templates in gee/ that differ from a fresh render (CI: must be empty)."""
    return render(check_only=True)


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
