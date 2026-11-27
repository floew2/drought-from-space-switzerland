"""Reading a drought from space: Swiss grassland in 2026, from Sentinel-2.

The package is laid out by responsibility:

    drought.core       configuration (validated), AOIs, grid, provenance, results
    drought.sources    swissEO (catalogue, image store, fetch), MeteoSwiss, parcels
    drought.mask       the PSF purity mask
    drought.products   series computed from the stored images
    drought.analysis   cleaning, season integral, cut detection, evaluation
    drought.viz        every figure
    drought.pipeline   the stages and the runner behind `drought run`

A short surface for notebooks:

    import drought as dr
    df = dr.load_site("linthebene")        cleaned per-parcel cut series
    dr.detect.detect_all(df)                any analysis module by its name
"""

import importlib as _importlib

# Short names for the modules a notebook uses most, imported on first use.
_MODULES = {
    "config": "drought.core.config", "provenance": "drought.core.provenance",
    "results": "drought.core.results", "aoi": "drought.core.aoi",
    "clean": "drought.analysis.clean", "detect": "drought.analysis.detect",
    "changepoint": "drought.analysis.changepoint", "evaluate": "drought.analysis.evaluate",
    "metrics": "drought.analysis.metrics", "labels": "drought.analysis.labels",
    "labelling": "drought.analysis.labelling", "io": "drought.analysis.io",
    "condition": "drought.analysis.condition", "mask": "drought.mask",
    "synthetic": "drought.synthetic", "demo": "drought.demo", "notebook": "drought.notebook",
    "cli": "drought.cli", "pipeline": "drought.pipeline",
}
_FROM_CONFIG = ("CFG", "available", "figures_dir", "write_provenance")

__all__ = [*_FROM_CONFIG, *_MODULES, "load_site"]


def __getattr__(name: str):
    """Import modules on first use.

    Lazy on purpose: some modules copy settings into constants when imported,
    so `import drought` must not import them before the `drought` command has
    applied its overrides.
    """
    if name in _MODULES:
        return _importlib.import_module(_MODULES[name])
    if name in _FROM_CONFIG:
        return getattr(_importlib.import_module("drought.core.config"), name)
    raise AttributeError(f"module 'drought' has no attribute {name!r}")


def load_site(site: str, verbose: bool = True, labelled: bool = False):
    """Read and clean one site's per-parcel cut series in one call."""
    from drought.analysis.clean import clean as _clean
    from drought.analysis.io import load_cuts
    return _clean(load_cuts(site, verbose=verbose, labelled=labelled), verbose=verbose)
