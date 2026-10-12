"""Configuration and path resolution.

One source of settings, loaded once. The previous notebook kept them as
module-level globals edited by hand, which meant a stale cell could leave a
parameter defined but unused — that happened three times, and once produced a
confidently wrong conclusion. Here the settings travel as an object, and
`provenance()` dumps them next to any result.
"""

from __future__ import annotations

import datetime as _dt
import json
import os as _os
import sys as _sys
import warnings as _warnings
from pathlib import Path
from types import SimpleNamespace

import yaml

_HERE = Path(__file__).resolve()
ROOT = _HERE.parents[2]          # <root>/src/drought/config.py -> <root>


class Config(SimpleNamespace):
    """Attribute access over the YAML, with nested dicts wrapped recursively."""

    def __init__(self, mapping: dict):
        super().__init__(**{
            k: Config(v) if isinstance(v, dict) else v
            for k, v in mapping.items()
        })

    def to_dict(self) -> dict:
        return {
            k: v.to_dict() if isinstance(v, Config) else v
            for k, v in vars(self).items()
        }

    def get(self, dotted: str, default=None):
        """cfg.get('detect.min_drop_abs')"""
        node = self
        for part in dotted.split("."):
            if not hasattr(node, part):
                return default
            node = getattr(node, part)
        return node


def load(path: str | Path | None = None) -> Config:
    path = Path(path) if path else ROOT / "config.yaml"
    with open(path) as f:
        return Config(yaml.safe_load(f))


CFG = load()
_REAL_PATHS = Config(CFG.paths.to_dict())
_PROFILE = "real"


def activate(profile: str = "real") -> None:
    """Switch every data and output path to a named profile.

    "real" uses `paths:`; "synthetic" uses `synthetic.paths:`, so synthetic
    inputs and outputs live in their own folders and can never be mistaken for,
    or overwrite, real ones. Settings other than paths are shared: the demo
    runs the same method with the same parameters.

    Args:
        profile: "real" or "synthetic".

    Raises:
        ValueError: For an unknown profile.
    """
    global _PROFILE
    if profile == "real":
        CFG.paths = Config(_REAL_PATHS.to_dict())
    elif profile == "synthetic":
        CFG.paths = Config(CFG.synthetic.paths.to_dict())
    else:
        raise ValueError(f"unknown profile {profile!r}; use 'real' or 'synthetic'")
    _PROFILE = profile


def active_profile() -> str:
    """The profile currently in force: "real" or "synthetic"."""
    return _PROFILE


# Modules that copy settings into module-level constants when first imported.
# An override applied after one of them is loaded would never reach it — the
# setting would look changed and change nothing, which is how a stale
# notebook cell once produced a confidently wrong conclusion. So it is refused.
_READ_AT_IMPORT = ("drought.mask.settings", "drought.trajectory", "drought.precip")


def _refuse_if_read(action: str) -> None:
    loaded = [m for m in _READ_AT_IMPORT if m in _sys.modules]
    if loaded:
        raise RuntimeError(
            f"cannot {action}: {', '.join(loaded)} already read the configuration "
            "when imported, so the change would not reach it. Change config.yaml, "
            "or apply overrides before importing the analysis modules (the "
            "`drought` command does this).")


def use_config(path: str | Path) -> None:
    """Replace the loaded configuration with another YAML file, in place.

    Args:
        path: The YAML file to load instead of config.yaml.

    Raises:
        RuntimeError: If an analysis module has already read the configuration.
    """
    global _REAL_PATHS, _PROFILE
    _refuse_if_read(f"load {path}")
    new = load(path)
    CFG.__dict__.clear()
    CFG.__dict__.update(vars(new))
    _REAL_PATHS = Config(CFG.paths.to_dict())
    _PROFILE = "real"


def override(dotted: str, value) -> None:
    """Set one existing setting, e.g. override("mask.purity_min", 0.8).

    Only existing keys can be set, so a typo fails instead of adding a setting
    nothing reads. Overriding the frozen cut-detector configuration is allowed
    but warned about: its validated F1 holds only for the frozen values.

    Args:
        dotted: Dotted key as in config.yaml.
        value: New value, already typed.

    Raises:
        KeyError: If the key does not exist.
        RuntimeError: If an analysis module has already read the configuration.
    """
    _refuse_if_read(f"override {dotted}")
    *parents, leaf = dotted.split(".")
    node = CFG
    for i, part in enumerate(parents):
        if not isinstance(getattr(node, part, None), Config):
            raise KeyError(f"no section {'.'.join(parents[:i + 1])!r} in the configuration")
        node = getattr(node, part)
    if not hasattr(node, leaf):
        known = ", ".join(sorted(vars(node)))
        raise KeyError(f"no setting {dotted!r}. Known in this section: {known}")
    if isinstance(getattr(node, leaf), Config):
        raise KeyError(f"{dotted!r} is a section, not a setting")
    if parents and parents[0] == "detect":
        _warnings.warn(f"{dotted} overridden: the cut detector configuration is frozen, "
                       "and its validated scores do not apply to other values",
                       stacklevel=2)
    setattr(node, leaf, value)


if _os.environ.get("DROUGHT_PROFILE"):
    activate(_os.environ["DROUGHT_PROFILE"])


# --------------------------------------------------------------------------- #
# Paths — site-aware, created on demand
# --------------------------------------------------------------------------- #
def _base(kind: str) -> Path:
    return ROOT / getattr(CFG.paths, kind)


def gee_dir(site: str) -> Path:
    """Where GEE exports for a site live. You drop the CSVs in here by hand."""
    p = _base("gee") / site
    p.mkdir(parents=True, exist_ok=True)
    return p


def figures_dir(site: str | None = None) -> Path:
    p = _base("figures") / site if site else _base("figures")
    p.mkdir(parents=True, exist_ok=True)
    return p


def tables_dir() -> Path:
    p = _base("tables")
    p.mkdir(parents=True, exist_ok=True)
    return p


def labels_path(site: str) -> Path:
    """One label file per site.

    Reference data is per-site on purpose: a detector validated on Linthebene
    has not been validated on Rheintal, and pooling the labels would hide that.
    """
    p = _base("labels")
    p.mkdir(parents=True, exist_ok=True)
    return p / f"cut_labels_{site}.json"


def parcels_dir() -> Path:
    """Yearly Nutzungsflaechen shapefiles. Searched recursively, since each
    year usually arrives in its own subfolder."""
    p = _base("parcels")
    p.mkdir(parents=True, exist_ok=True)
    return p


def station_dir() -> Path:
    """Cached MeteoSwiss daily series (or their synthetic stand-ins)."""
    p = _base("station")
    p.mkdir(parents=True, exist_ok=True)
    return p


def mask_path(name: str = "grassland_mask_sg.tif") -> Path:
    p = _base("masks")
    p.mkdir(parents=True, exist_ok=True)
    return p / name


def site_cfg(site: str) -> Config:
    if not hasattr(CFG.sites, site):
        known = ", ".join(vars(CFG.sites))
        raise KeyError(f"unknown site {site!r}. Known sites: {known}")
    return getattr(CFG.sites, site)


def available(site: str) -> dict:
    """What data actually exists for a site, so a notebook can say so plainly
    instead of failing several cells later."""
    d = gee_dir(site)
    cuts = sorted(d.glob("cuts_*.csv"))
    # The trajectory export covers all sites in one file, so it may sit in the
    # site folder or in data/gee/ — look in both.
    traj = sorted(d.glob("*trajectory*.csv")) + \
           sorted((_base("gee")).glob("*trajectory*.csv"))
    return {
        "site": site,
        "gee_dir": str(d.relative_to(ROOT)),
        "cut_years": sorted(int(f.stem.split("_")[-1]) for f in cuts
                            if f.stem.split("_")[-1].isdigit()),
        "trajectory_files": [f.name for f in traj],
        "labels": labels_path(site).exists(),
    }


def provenance(extra: dict | None = None) -> dict:
    """Settings + timestamp, to be written beside any result.

    A number without its configuration is not reproducible, and this one has
    changed six times.
    """
    out = {
        "generated": _dt.datetime.now().isoformat(timespec="seconds"),
        "config": CFG.to_dict(),
    }
    if extra:
        out.update(extra)
    return out


def write_provenance(name: str, extra: dict | None = None) -> Path:
    p = tables_dir() / f"{name}_provenance.json"
    with open(p, "w") as f:
        json.dump(provenance(extra), f, indent=1, default=str)
    return p
