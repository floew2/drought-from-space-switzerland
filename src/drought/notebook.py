"""Run notebooks/walkthrough.ipynb on synthetic data and commit-proof its outputs.

    drought notebook            # execute in place (the committed version)
    drought notebook --check    # execute a copy; fail if any cell errors (CI)

The committed notebook carries the text outputs of a SYNTHETIC run (numbers
and tables), so GitHub shows what each step returns. Figures are dropped:
synthetic figures are never committed, and they appear when you run it.
Before saving, local paths are replaced by repository-relative ones. `tests/test_notebook.py` refuses a committed notebook that was run on
real data, or that leaks a local path.
"""

from __future__ import annotations

import re
import tempfile
from pathlib import Path

from .config import ROOT

NOTEBOOK = ROOT / "notebooks" / "walkthrough.ipynb"
_HOME = re.compile(r"(/Users/|/home/|C:\\\\Users\\\\)[^\s'\"]*?(?=data/|output/|notebooks/)")


def _scrub(text: str) -> str:
    text = text.replace(str(ROOT) + "/", "")
    return _HOME.sub("", text)


def _clean(nb) -> None:
    for c in nb.cells:
        if c.cell_type != "code":
            continue
        c.metadata.pop("execution", None)
        c.outputs = [o for o in c.get("outputs", [])
                     if "image/png" not in o.get("data", {})]
        for o in c.outputs:
            if o.output_type == "stream":
                o.text = _scrub(o.text)
            elif "data" in o and "text/plain" in o.data:
                o.data["text/plain"] = _scrub(o.data["text/plain"])


def run(check: bool = False, path: Path = NOTEBOOK) -> Path:
    """Execute the walkthrough (synthetic profile) and save it, cleaned."""
    try:
        import nbformat
        from nbclient import NotebookClient
    except ImportError as e:
        raise SystemExit('needs the notebook extra: pip install -e ".[notebook]"') from e
    nb = nbformat.read(path, as_version=4)
    src = next(c.source for c in nb.cells if c.cell_type == "code")
    if not re.search(r"^REAL\s*=\s*False", src, re.M):
        raise SystemExit("walkthrough: set REAL = False before running it for the repository")
    NotebookClient(nb, timeout=900, kernel_name="python3",
                   resources={"metadata": {"path": str(path.parent)}}).execute()
    _clean(nb)
    out = Path(tempfile.mkdtemp()) / path.name if check else path
    nbformat.write(nb, out)
    print(f"  {path.relative_to(ROOT)}: all cells ran" + ("" if check else "; saved"))
    return out
