"""Drought 2026 — Swiss grassland and water analysis from Sentinel-2.

Import the pieces you need:

    from drought import config, io, clean, detect, metrics, evaluate
    from drought import figures, labels, labelling, precip

Or the convenience surface:

    import drought as dr
    df = dr.load_site("linthebene")
"""

import importlib as _importlib

_SUBMODULES = ("changepoint", "clean", "cli", "config", "demo", "detect", "evaluate",
               "figures", "gallery", "gee", "io", "labelling", "labels", "mask",
               "mask_figures", "metrics", "precip", "provenance", "results", "synthetic",
               "trajectory")
# config.provenance() (settings + timestamp) is reached as drought.config.provenance;
# drought.provenance is the module that records and verifies files.
_FROM_CONFIG = ("CFG", "available", "figures_dir", "write_provenance")

__all__ = [*_FROM_CONFIG, *_SUBMODULES, "load_site"]


def __getattr__(name: str):
    """Import submodules on first use.

    Lazy on purpose: some modules copy settings into constants when imported,
    so `import drought` must not import them before the `drought` command has
    applied its overrides.
    """
    if name in _SUBMODULES:
        return _importlib.import_module(f".{name}", __name__)
    if name in _FROM_CONFIG:
        return getattr(_importlib.import_module(".config", __name__), name)
    raise AttributeError(f"module 'drought' has no attribute {name!r}")


def load_site(site: str, verbose: bool = True):
    """Read and clean one site's per-parcel cut data in one call."""
    from .clean import clean as _clean
    from .io import load_cuts
    return _clean(load_cuts(site, verbose=verbose), verbose=verbose)
