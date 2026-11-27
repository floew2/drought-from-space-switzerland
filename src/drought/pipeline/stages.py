"""The stages of the pipeline: what each needs, what it makes, when it is done.

    fetch        swissEO mosaics (and VHI maps) for an AOI -> data/scenes/
    mask         purity mask and parcel ids from the parcel data -> data/masks/, data/parcels/
    condition    site-mean index series                  -> data/series/<site>/condition.csv
    sensitivity  anomaly per purity class                -> sensitivity_<regime>.csv
    anomaly      per-field z-scores                       -> anomaly.csv
    cuts         per-parcel series for cut detection     -> cuts_<year>[_labelled].csv
    vhi          swisstopo's VHI on the same cells       -> vhi_series.csv, vhi_parcels.csv
    evaluate     cut detector scored against the labels  -> results.json
    precip       MeteoSwiss station series               -> data/station/
    report       every figure and headline number        -> output/

A stage is done when its outputs exist and their provenance records match the
current settings and mask (core/provenance.py), so a changed setting re-runs
exactly the stages it affects.
"""

from __future__ import annotations

import logging
from collections.abc import Callable
from dataclasses import dataclass, field

from drought.core import provenance
from drought.core.config import CFG, labels_path, mask_path, station_dir, years

__all__ = ["STAGES", "Stage", "sites_for"]

log = logging.getLogger(__name__)


@dataclass(frozen=True)
class Stage:
    name: str
    help: str
    needs: tuple[str, ...]
    per_site: bool
    run: Callable[..., object]
    done: Callable[..., bool]
    default_sites: Callable[[], list[str]] = field(default=lambda: [])


# --------------------------------------------------------------------------- #
# Which sites a stage runs for by default
# --------------------------------------------------------------------------- #
def _all_sites() -> list[str]:
    return list(dict.fromkeys([*CFG.condition.sites, *CFG.cuts.sites, *CFG.anomaly.sites]))


def _labelled_sites() -> list[str]:
    return [s for s in CFG.cuts.sites if labels_path(s).exists()]


def sites_for(stage: Stage, requested: list[str] | None) -> list[str]:
    """The sites a stage runs for: those requested, else its defaults."""
    return list(requested) if requested else stage.default_sites()


# --------------------------------------------------------------------------- #
# The stages
# --------------------------------------------------------------------------- #
def _fetch(site: str, years_: list[int] | None = None, workers: int | None = None, **_):
    from drought.sources import fetch
    failed = fetch.run([site], years_ or years(), workers=workers)
    if failed:
        raise RuntimeError(f"{site}: {len(failed)} mosaic(s) not downloaded; run again")


def _fetch_done(site: str, years_: list[int] | None = None, **_) -> bool:
    from drought.sources import store
    ys = years_ or years()
    ok = all(store.complete(site, y) for y in ys)
    if int(CFG.project.current_year) in ys:
        ok &= store.complete(site, int(CFG.project.current_year), "vhi")
    return ok


def _mask(_site=None, **_):
    from drought.core import results
    from drought.mask import build_mask, psf_params
    from drought.sources import parcels
    master = build_mask()
    results.update("mask", {"funnel": master.attrs["funnel"], "response": psf_params()},
                   inputs=[mask_path()])
    parcels.write()


def _mask_done(_site=None, **_) -> bool:
    from drought.sources import parcels
    return provenance.is_current(mask_path(), "mask") and parcels.path().exists()


def _product(module: str):
    def run(site: str, **_):
        import importlib
        importlib.import_module(f"drought.products.{module}").run(site)
    return run


def _condition_done(site: str, **_) -> bool:
    from drought.products import condition
    return provenance.is_current(condition.output(site), "condition")


def _sensitivity_done(site: str, **_) -> bool:
    from drought.products import sensitivity
    return all(provenance.is_current(sensitivity.output(site, r), "sensitivity")
               for r in CFG.analysis.regimes)


def _anomaly_done(site: str, **_) -> bool:
    from drought.products import anomaly
    return provenance.is_current(anomaly.output(site), "anomaly")


def _cuts_done(site: str, **_) -> bool:
    from drought.products import cuts
    lab = cuts.labelled_ids(site)
    return all(provenance.is_current(cuts.output(site, y), "cuts") and
               (y not in lab or provenance.is_current(cuts.output(site, y, True), "cuts"))
               for y in years())


def _vhi_done(site: str, **_) -> bool:
    from drought.products import vhi
    return all(provenance.is_current(vhi.output(site, w), "vhi") for w in ("series", "parcels"))


def _evaluate(site: str, **_):
    from drought.pipeline.report import evaluate
    evaluate(site)


def _precip(_site=None, **_):
    from drought.sources import meteoswiss
    meteoswiss.download(station_dir(), CFG.precip.station_name)


def _precip_done(_site=None, **_) -> bool:
    files = sorted(station_dir().glob("ogd-smn_*_d_*.csv"))
    return bool(files) and all(provenance.is_current(f, "station") for f in files)


def _report(_site=None, sites=None, **_):
    from drought.pipeline.report import report
    report(sites)


def _never(*_a, **_k) -> bool:
    return False


STAGES: dict[str, Stage] = {s.name: s for s in [
    Stage("fetch", "download swissEO for an AOI (no mask needed)", (), True,
          _fetch, _fetch_done, _all_sites),
    Stage("mask", "purity mask and stable parcel ids from the parcel data", (), False,
          _mask, _mask_done),
    Stage("condition", "site-mean index series per management system", ("fetch", "mask"),
          True, _product("condition"), _condition_done, lambda: list(CFG.condition.sites)),
    Stage("sensitivity", "anomaly per purity class", ("fetch", "mask"), True,
          _product("sensitivity"), _sensitivity_done, lambda: list(CFG.condition.sites)),
    Stage("anomaly", "per-field z-scores against each field's own baseline",
          ("fetch", "mask"), True, _product("anomaly"), _anomaly_done,
          lambda: list(CFG.anomaly.sites)),
    Stage("cuts", "per-parcel series for cut detection", ("fetch", "mask"), True,
          _product("cuts"), _cuts_done, lambda: list(CFG.cuts.sites)),
    Stage("vhi", "swisstopo's Vegetation Health Index on the same cells", ("fetch", "mask"),
          True, _product("vhi"), _vhi_done, lambda: list(CFG.condition.sites)),
    Stage("evaluate", "score the cut detector against the hand labels", ("cuts",), True,
          _evaluate, _never, _labelled_sites),
    Stage("precip", "MeteoSwiss station series", (), False, _precip, _precip_done),
    Stage("report", "every figure and headline number",
          ("condition", "sensitivity", "anomaly", "vhi", "evaluate", "precip"), False,
          _report, _never),
]}
