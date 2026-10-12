"""Drought 2026 — Swiss grassland and water analysis from Sentinel-2.

Import the pieces you need:

    from drought import config, io, clean, detect, metrics, evaluate
    from drought import figures, labels, labelling, precip

Or the convenience surface:

    import drought as dr
    df = dr.load_site("linthebene")
"""

from . import (
    changepoint,
    clean,
    config,
    detect,
    evaluate,
    figures,
    io,
    labelling,
    labels,
    mask,
    mask_figures,
    metrics,
    precip,
    trajectory,
)
from .config import CFG, available, figures_dir, provenance, write_provenance

__all__ = ["CFG", "available", "figures_dir", "provenance", "write_provenance",
           "changepoint", "clean", "config", "detect", "evaluate", "figures",
           "io", "labels", "labelling", "mask", "mask_figures", "metrics", "precip",
           "trajectory",
           "load_site"]


def load_site(site: str, verbose: bool = True):
    """Read and clean one site's per-parcel cut data in one call."""
    from .clean import clean as _clean
    from .io import load_cuts
    return _clean(load_cuts(site, verbose=verbose), verbose=verbose)
