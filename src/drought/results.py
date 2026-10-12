"""output/results.json: every headline number, in one artifact, with its provenance.

The README, the deck and METHOD cite this file rather than numbers copied from a
terminal. Each step that produces a headline number writes its own section:

    mask            master-set funnel, per-category totals, sensor response
    mask_area       observed / polygon area at each purity threshold
    condition       season NDMI integral per site and management system
    precipitation   cumulative P and P - ET0 against the 1991-2020 median
    detector.<site> precision, recall, F1, per-year recall, changepoint F1

Every section records the settings hash and mask it was computed under. A
section computed under other settings than the current ones is stale and is
removed when the file is next written, with a notice, so the file never mixes
vintages. `load()` refuses a file containing any stale section.

Under the synthetic profile the file is output/synthetic/results.json.
"""

from __future__ import annotations

import datetime as dt
import json
from pathlib import Path

from .config import CFG, ROOT, active_profile

__all__ = ["load", "path", "update"]


def path() -> Path:
    """results.json for the active profile: beside the tables folder."""
    return (ROOT / CFG.paths.tables).parent / "results.json"


def _stamp() -> dict:
    from .provenance import mask_fingerprint, settings, settings_hash

    fp = mask_fingerprint()
    return {"settings_hash": settings_hash(settings("results")),
            "mask_sha256": fp["sha256"] if fp else None}


def update(section: str, data: dict, inputs: list | tuple = (), merge: bool = False) -> Path:
    """Write one section, stamped with the current settings and mask.

    Args:
        section: Section name, e.g. "condition" or "detector.linthebene".
        data: The numbers. Must be JSON-serialisable.
        inputs: Files the numbers were computed from. Each one's provenance
            record is referenced by its hash.
        merge: Add `data`'s keys to the existing section instead of replacing
            it (one site at a time, for example).

    Returns:
        Path of results.json.
    """
    from .provenance import record, sidecar_path

    p = path()
    p.parent.mkdir(parents=True, exist_ok=True)
    doc = json.loads(p.read_text()) if p.exists() else {"sections": {}}
    now = _stamp()
    for name in list(doc["sections"]):
        old = doc["sections"][name]
        if {k: old.get(k) for k in now} != now:
            print(f"  results.json: removed stale section {name!r} "
                  f"(computed under settings {old.get('settings_hash')}, "
                  f"mask {old.get('mask_sha256')})")
            del doc["sections"][name]
    refs = []
    for f in inputs:
        side = sidecar_path(f)
        rec = json.loads(side.read_text()) if side.exists() else {}
        refs.append({"file": rec.get("file", str(f)), "file_sha256": rec.get("file_sha256"),
                     "recorded": rec.get("created")})
    if merge and section in doc["sections"]:
        old = doc["sections"][section]
        data = {**old["data"], **data}
        refs = [r for r in old["inputs"] if r not in refs] + refs
    doc["sections"][section] = {**now, "data": data, "inputs": refs,
                                "computed": dt.datetime.now().astimezone()
                                .isoformat(timespec="seconds")}
    doc["profile"] = active_profile()
    doc["provenance"] = record("results")
    doc["sections"] = dict(sorted(doc["sections"].items()))
    p.write_text(json.dumps(doc, indent=1, default=_jsonable))
    return p


def _jsonable(o):
    import numpy as np
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, (tuple, set)):
        return list(o)
    return str(o)


def load(strict: bool = True) -> dict:
    """Read results.json, refusing it if any section is stale.

    Args:
        strict: Raise on a stale section (default) rather than returning it.

    Returns:
        The document.

    Raises:
        FileNotFoundError: If no results have been written.
        RuntimeError: If a section was computed under other settings or another mask.
    """
    doc = json.loads(path().read_text())
    now = _stamp()
    stale = [n for n, s in doc["sections"].items() if {k: s.get(k) for k in now} != now]
    if stale and strict:
        raise RuntimeError(f"results.json: sections {stale} were computed under other "
                           "settings or another mask; re-run them")
    return doc
