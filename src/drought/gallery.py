"""Render the model figure for docs/gallery.md.

    python -m drought.gallery

docs/figures/model/   figures that use no data at all (the sensor response)
docs/figures/real/    results from the real data, copied in by hand from a
                      real run, never generated here (see the README there)

Synthetic figures are never committed: `make demo` writes them to
output/synthetic/figures/, for anyone who wants to see what each step does.
"""

from __future__ import annotations

from pathlib import Path

from .config import ROOT

DOCS = ROOT / "docs" / "figures"
MAX_PX = 1200


def _shrink(src: Path, dst: Path) -> None:
    from PIL import Image

    im = Image.open(src)
    im.thumbnail((MAX_PX, MAX_PX))
    dst.parent.mkdir(parents=True, exist_ok=True)
    # No palette quantisation: it shifts line colours, and the legends key on them.
    im.convert("RGB").save(dst, optimize=True)


def build() -> list[Path]:
    """Render the model figure. Returns the written paths."""
    import matplotlib

    matplotlib.use("Agg")
    from . import mask_figures

    out = []
    model = DOCS / "model"
    model.mkdir(parents=True, exist_ok=True)
    tmp = model / "_psf_full.png"
    mask_figures.figure_psf(out=tmp)
    _shrink(tmp, model / "psf.png")
    tmp.unlink()
    out.append(model / "psf.png")

    for p in out:
        print(f"  {p.relative_to(ROOT)}  ({p.stat().st_size / 1e3:.0f} KB)")
    return out


if __name__ == "__main__":
    build()
