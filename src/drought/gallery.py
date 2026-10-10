"""Build the committed figures for docs/gallery.md.

    python -m drought.gallery        # after `make demo`

Two sources, kept apart on disk and in the gallery:

docs/figures/model/   figures that use no data at all (the sensor response),
                      rendered here without the synthetic stamp
docs/figures/demo/    the demo's synthetic figures, downscaled; they keep the
                      SYNTHETIC stamp they were saved with
docs/figures/real/    results from the real data — copied in by hand from a
                      real run, never generated here (see the README there)
"""

from __future__ import annotations

from pathlib import Path

from .config import ROOT

DOCS = ROOT / "docs" / "figures"
DEMO_PICKS = [
    "linthebene/vi_trajectory_NDMI_intensive.png",
    "rheintal/vi_trajectory_NDMI_extensive.png",
    "station/precip_syn_2026.png",
    "mask/parcel_size.png",
    "mask/area_check_all.png",
    "mask/purity_explainer_linthebene_intensive.png",
    "mask/purity_map_linthebene_intensive.png",
    "linthebene/last_cut_intensive.png",
]
MAX_PX = 1200


def _shrink(src: Path, dst: Path) -> None:
    from PIL import Image

    im = Image.open(src)
    im.thumbnail((MAX_PX, MAX_PX))
    dst.parent.mkdir(parents=True, exist_ok=True)
    # No palette quantisation: it shifts line colours, and the legends key on them.
    im.convert("RGB").save(dst, optimize=True)


def build() -> list[Path]:
    """Render the model figure and copy the demo picks. Returns written paths."""
    import matplotlib

    matplotlib.use("Agg")
    from . import mask
    from .config import CFG

    out = []
    model = DOCS / "model"
    model.mkdir(parents=True, exist_ok=True)
    tmp = model / "_psf_full.png"
    mask.figure_psf(out=tmp)
    _shrink(tmp, model / "psf.png")
    tmp.unlink()
    out.append(model / "psf.png")

    demo_root = ROOT / CFG.synthetic.paths.figures
    missing = [p for p in DEMO_PICKS if not (demo_root / p).exists()]
    if missing:
        raise SystemExit(f"run `make demo` first; missing: {missing}")
    for p in DEMO_PICKS:
        dst = DOCS / "demo" / p.replace("/", "_")
        _shrink(demo_root / p, dst)
        out.append(dst)
    for p in out:
        print(f"  {p.relative_to(ROOT)}  ({p.stat().st_size / 1e3:.0f} KB)")
    return out


if __name__ == "__main__":
    build()
