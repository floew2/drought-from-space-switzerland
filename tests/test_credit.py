"""Every saved figure names its data sources and the project."""

import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from drought.core.config import CFG  # noqa: E402
from drought.viz import credit  # noqa: E402

SRC = Path(__file__).resolve().parents[1] / "src" / "drought"


def test_no_figure_is_saved_around_the_credit():
    direct = [f"{p.relative_to(SRC)}:{i}"
              for p in SRC.rglob("*.py") if p.name != "credit.py" and "synthetic" not in p.parts
              for i, ln in enumerate(p.read_text().splitlines(), 1)
              if re.search(r"\.savefig\(", ln)]
    assert not direct, f"save through drought.viz.credit.save: {direct}"


def test_the_line_carries_the_required_wording():
    s2 = credit.line("s2", "parcels")
    assert "© swisstopo" in s2 and "contains modified Copernicus Sentinel data" in s2
    assert "Kanton St.Gallen" in s2 and s2.endswith(CFG.project.credit)
    assert "MeteoSwiss" in credit.line("meteo")


def test_save_stamps_once(tmp_path):
    fig, _ = plt.subplots()
    credit.save(fig, tmp_path / "a.png", "meteo")
    credit.save(fig, tmp_path / "b.png", "meteo")
    assert sum(t.get_text() == credit.line("meteo") for t in fig.texts) == 1
    assert (tmp_path / "b.png").stat().st_size > 0
    plt.close(fig)
