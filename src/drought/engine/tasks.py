"""Start table exports to Google Drive, wait, download and register what they wrote.

Export names carry the CONFIG_HASH of the settings and mask they ran under
(`..._cfg<hash>`), exactly as the rendered Code Editor scripts did, so
`provenance.register` can refuse an export from other settings. Registering
happens as part of the download: a file arrives with its record or not at all.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path

from ..config import CFG, ROOT
from . import session

__all__ = ["Export", "config_hash", "fetch", "run", "start", "wait"]

POLL_S = 30


def config_hash() -> str:
    """The hash of every Earth Engine setting plus the local mask."""
    from ..gee import _remember, expected_hash

    h = expected_hash(with_mask=True)
    _remember(ROOT / "output" / "gee", h)      # so a mismatch can name what changed
    return h


@dataclass
class Export:
    """One table export: what to compute, where it lands, how it is checked."""
    collection: object            # ee.FeatureCollection
    prefix: str                   # file name without .csv, including _cfg<hash>
    dest_dir: Path                # local folder for the CSV
    kind: str | None              # provenance kind; None = not registered (smoke runs)
    selectors: list[str] | None = None
    task: object = None

    @property
    def filename(self) -> str:
        return f"{self.prefix}.csv"

    @property
    def local(self) -> Path:
        return self.dest_dir / f"{self.prefix}.csv"


def start(exp: Export) -> Export:
    """Submit one export to Google Drive (folder gee.drive_folder)."""
    ee = session.ee()
    exp.task = ee.batch.Export.table.toDrive(
        collection=exp.collection, description=exp.prefix[:100],
        folder=str(CFG.gee.drive_folder), fileNamePrefix=exp.prefix, fileFormat="CSV",
        selectors=exp.selectors)
    exp.task.start()
    print(f"  started {exp.prefix}  (task {exp.task.id})")
    return exp


def wait(exports: list[Export], poll_s: int = POLL_S) -> None:
    """Block until every task has finished; raise if any failed.

    Prints each state change and, when Earth Engine reports it, the EECU cost.
    """
    last = {}
    pending = [e for e in exports if e.task is not None]
    while pending:
        for e in list(pending):
            st = e.task.status()
            state = st.get("state")
            if last.get(e.prefix) != state:
                eecu = st.get("batch_eecu_usage_seconds")
                extra = f", {eecu:,.0f} EECU-s" if eecu else ""
                print(f"  {e.prefix}: {state}{extra}")
                last[e.prefix] = state
            if state == "COMPLETED":
                pending.remove(e)
            elif state in ("FAILED", "CANCELLED", "CANCEL_REQUESTED"):
                raise RuntimeError(f"{e.prefix}: {state} — {st.get('error_message', '')}")
        if pending:
            time.sleep(poll_s)


def fetch(exp: Export) -> Path:
    """Download one finished export and register its provenance."""
    from ..provenance import register

    found = session.drive_find(exp.filename)
    if not found:
        raise FileNotFoundError(f"{exp.filename} not found in Google Drive")
    session.drive_download(found[0]["id"], exp.local)     # newest, if Drive holds two
    if exp.kind:
        register(exp.local, exp.kind)
    print(f"  -> {exp.local.relative_to(ROOT) if exp.local.is_relative_to(ROOT) else exp.local}"
          + ("  (registered)" if exp.kind else "  (smoke run: not registered)"))
    return exp.local


def run(exports: list[Export], dry_run: bool = False) -> list[Path]:
    """Check the mask asset, then start, wait for, download and register each export.

    Args:
        exports: The exports to run.
        dry_run: Print what would run and stop.

    Returns:
        Local paths of the downloaded CSVs.
    """
    from .imagery import check_mask_asset

    sha = check_mask_asset()
    print(f"  mask asset {CFG.gee.assets.mask} matches the local mask ({sha})")
    if dry_run:
        for e in exports:
            print(f"  would export {e.prefix}  ->  {e.dest_dir}")
        return []
    todo = []
    for e in exports:
        if e.local.exists():
            print(f"  {e.local.name} already here; skipped")
            continue
        if session.drive_find(e.filename):
            print(f"  {e.prefix} already exported; downloading")
        else:
            start(e)
        todo.append(e)
    wait(todo)
    return [fetch(e) for e in todo]
