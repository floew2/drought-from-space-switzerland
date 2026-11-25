"""Run a stage, or a stage and everything it needs, skipping what is done.

    runner.run("cuts", sites=["linthebene"])               cuts and what it needs
    runner.run("cuts", sites=["linthebene"], deps=False)   cuts alone
    runner.plan("report")                                  what would run, in order

A stage is skipped when it is done (its outputs exist and match the current
settings and mask); `force=True` re-runs the requested stage regardless.
"""

from __future__ import annotations

import logging

from drought.pipeline.stages import STAGES, Stage, sites_for

__all__ = ["order", "plan", "run"]

log = logging.getLogger(__name__)


def order(target: str) -> list[Stage]:
    """The target and everything it needs, dependencies first."""
    if target not in STAGES:
        raise KeyError(f"unknown stage {target!r}; stages: {', '.join(STAGES)}")
    seen: list[str] = []

    def visit(name: str, path: tuple[str, ...] = ()) -> None:
        if name in path:
            raise RuntimeError(f"stage cycle: {' -> '.join((*path, name))}")
        for dep in STAGES[name].needs:
            visit(dep, (*path, name))
        if name not in seen:
            seen.append(name)

    visit(target)
    return [STAGES[n] for n in seen]


def _jobs(stage: Stage, sites: list[str] | None) -> list[str | None]:
    jobs: list[str | None] = list(sites_for(stage, sites)) if stage.per_site else [None]
    return jobs


def plan(target: str, sites: list[str] | None = None, deps: bool = True,
         force: bool = False, **kw) -> list[tuple[str, str | None, bool]]:
    """(stage, site, will_run) for every job, in the order they would run."""
    stages = order(target) if deps else [STAGES[target]]
    out = []
    for st in stages:
        for site in _jobs(st, sites):
            will = (force and st.name == target) or not st.done(site, sites=sites, **kw)
            out.append((st.name, site, will))
    return out


def run(target: str, sites: list[str] | None = None, deps: bool = True,
        force: bool = False, **kw) -> list[tuple[str, str | None]]:
    """Run `target` (and, with `deps`, whatever it needs that is not done).

    Args:
        target: A stage name (pipeline.stages.STAGES).
        sites: AOIs to run per-site stages for (default: each stage's own).
        deps: Also run the stages it needs, if they are not done.
        force: Re-run `target` even if it is done.
        **kw: Passed to the stages: years_ and workers for fetch.

    Returns:
        The (stage, site) jobs that ran.
    """
    ran = []
    for name, site, will in plan(target, sites, deps, force, **kw):
        label = f"{name} {site}" if site else name
        if not will:
            log.info("== %s: up to date", label)
            continue
        log.info("== %s %s", label, "=" * max(1, 60 - len(label)))
        STAGES[name].run(site, sites=sites, **kw)
        ran.append((name, site))
    return ran
