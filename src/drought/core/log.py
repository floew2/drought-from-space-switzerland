"""Logging: the package logs through `logging`; the command line prints it plainly.

Library modules call `logging.getLogger(__name__)` and never configure
anything. `setup` is called once by the `drought` command (and the demo), and
sends the package's messages to stdout without decoration, so the output
reads like a report. In a notebook, call `setup()` to see the same.
"""

from __future__ import annotations

import contextlib
import logging
import sys
from collections.abc import Iterator

__all__ = ["quiet", "setup"]


def setup(level: int = logging.INFO) -> None:
    """Route the `drought` loggers to stdout, message only. Idempotent."""
    root = logging.getLogger("drought")
    root.setLevel(level)
    if not any(getattr(h, "_drought", False) for h in root.handlers):
        h = logging.StreamHandler(sys.stdout)
        h.setFormatter(logging.Formatter("%(message)s"))
        h._drought = True                            # type: ignore[attr-defined]
        root.addHandler(h)
        root.propagate = False


@contextlib.contextmanager
def quiet(level: int = logging.WARNING) -> Iterator[None]:
    """Hold the package's messages below `level` for the duration."""
    root = logging.getLogger("drought")
    old = root.level
    root.setLevel(max(level, old))
    try:
        yield
    finally:
        root.setLevel(old)
