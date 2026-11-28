"""Every saved figure carries its data sources and the project credit.

    credit.save(fig, path, "s2", "parcels")

writes one grey line below the figure, outside the plot area, so it can never
collide with tick labels or captions, and saves with a tight bounding box.
The wording of each source is the attribution its licence asks for.
"""

from __future__ import annotations

from pathlib import Path

from drought.core.config import CFG, years

__all__ = ["SOURCES", "line", "save"]


def _s2() -> str:
    ys = years()
    span = f"{ys[0]}–{ys[-1]}" if len(ys) > 1 else str(ys[0])
    return f"© swisstopo (swissEO S2-SR), contains modified Copernicus Sentinel data {span}"


SOURCES = {
    "s2": _s2,
    "vhi": lambda: "VHI © swisstopo, MeteoSwiss (swissEO)",
    "parcels": lambda: "parcels © Kanton St.Gallen (Nutzungsflächen)",
    "meteo": lambda: "Source: MeteoSwiss (Open Data, CC BY 4.0)",
}


def line(*sources: str) -> str:
    """The credit line: the named sources, then the project credit."""
    return " · ".join([*(SOURCES[s]() for s in sources), str(CFG.project.credit)])


def save(fig, path: str | Path, *sources: str, dpi: int = 180) -> Path:
    """Add the credit line below the figure and save it.

    Args:
        fig: The matplotlib figure.
        path: Where to write it.
        *sources: Keys of SOURCES for the data the figure shows.
        dpi: Resolution.

    Raises:
        KeyError: For an unknown source.
    """
    text = line(*sources)
    if not getattr(fig, "_credit", None):
        fig._credit = fig.text(0.005, -0.005, text, transform=fig.transFigure,
                               ha="left", va="top", fontsize=6.5, color="#8b949e")
    fig.savefig(path, dpi=dpi, facecolor=fig.get_facecolor(), bbox_inches="tight",
                pad_inches=0.12)
    return Path(path)
