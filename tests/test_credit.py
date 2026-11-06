"""Every saved figure carries the project credit (config: project.credit).

Two checks. A static one over the whole package: every function that calls
`savefig` calls `credit(fig)` first, so a new figure function cannot be added
without it. And a runtime one: saving real figures through the public entry
points leaves the credit text on the canvas, exactly once.
"""

from __future__ import annotations

import ast
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

from drought import config  # noqa: E402

PKG = Path(config.__file__).parent
# Monkeypatches savefig itself (the synthetic stamp), so it is not a figure.
NOT_FIGURE_CODE = {PKG / "synthetic" / "__init__.py"}


def _calls(node, name):
    for n in ast.walk(node):
        if isinstance(n, ast.Call):
            f = n.func
            if (isinstance(f, ast.Attribute) and f.attr == name) or \
               (isinstance(f, ast.Name) and f.id == name):
                yield n


def _figure_functions():
    for path in sorted(PKG.rglob("*.py")):
        if path in NOT_FIGURE_CODE:
            continue
        tree = ast.parse(path.read_text())
        for fn in ast.walk(tree):
            if isinstance(fn, (ast.FunctionDef, ast.AsyncFunctionDef)):
                saves = list(_calls(fn, "savefig"))
                if saves:
                    yield path, fn, saves


def test_every_savefig_is_preceded_by_credit():
    found = list(_figure_functions())
    assert found, "no figure functions found; the scan is broken"
    missing = []
    for path, fn, saves in found:
        credits = [c.lineno for c in _calls(fn, "credit")]
        for s in saves:
            if not any(c < s.lineno for c in credits):
                missing.append(f"{path.relative_to(PKG)}:{s.lineno} in {fn.name}()")
    assert not missing, "savefig without credit(fig) before it:\n  " + "\n  ".join(missing)


def test_credit_text_and_position_follow_config():
    fig = plt.figure()
    config.credit(fig)
    config.credit(fig)                       # idempotent
    texts = [t for t in fig.texts if t.get_text() == config.CFG.project.credit]
    assert len(texts) == 1
    t = texts[0]
    assert t.get_position() == (0.995, 0.006)
    assert (t.get_ha(), t.get_va()) == ("right", "bottom")
    assert t.get_fontsize() == 7
    assert t.get_color() == "#8a8a8a" and t.get_alpha() == 0.85
    plt.close(fig)


@pytest.fixture
def saved(monkeypatch):
    """Record, for every savefig, whether the credit is on the canvas."""
    seen = []
    original = Figure.savefig

    def spy(self, *args, **kwargs):
        seen.append(sum(t.get_text() == config.CFG.project.credit for t in self.texts))
        return original(self, *args, **kwargs)

    monkeypatch.setattr(Figure, "savefig", spy)
    return seen


def test_entry_points_stamp_credit(tmp_path, saved):
    from drought import mask_figures, precip

    mask_figures.figure_psf(out=tmp_path / "psf.png")
    config.activate("synthetic")
    try:
        daily, abbr, lat = precip.load_local(config.station_dir())
        precip.main(daily=daily, abbr=abbr, lat=lat, station_label="Synthetic station")
    finally:
        config.activate("real")
    plt.close("all")
    assert len(saved) >= 2
    assert saved == [1] * len(saved), f"credit count per saved figure: {saved}"
