"""Mask settings, read once from config.yaml.

Everything the mask modules need from the configuration lives here, so a
setting has one name and one source. Edit config.yaml, not this file.
"""

from __future__ import annotations

from ..config import CFG

PATTERN = "nutzungen_*.shp"          # yearly deliveries, searched recursively
CRS = "EPSG:2056"                    # CH1903+ / LV95; asserted on load

RES = CFG.mask.res                   # 10 m, the Sentinel-2 grid
SUBPIXEL = CFG.mask.subpixel         # sub-cells per side: 4 -> 2.5 m
# Polygons below this are cadastral offcuts that would corrupt the overlap
# matching. NOT an observability filter: purity decides that (METHOD §2.1).
SLIVER_FLOOR_M2 = CFG.mask.sliver_floor_m2
PSF_MTF = CFG.mask.psf_mtf_nyquist
PSF_MODEL = getattr(CFG.mask, "psf_model", "box_gauss")
PURITY_MIN = CFG.mask.purity_min
# Share of a field that must lie inside a year's category to count as held.
MIN_OVERLAP = getattr(CFG.mask, "min_overlap", 0.8)
# Years the category must be held in; 0 = all.
MIN_YEARS_SAME = getattr(CFG.mask, "min_years_same", 0)
# Simplify before the overlay. Measured and off: it costs more than it saves.
SIMPLIFY_M = getattr(CFG.mask, "simplify_m", 0.0)

_clip = getattr(CFG.mask, "clip_bbox_lv95", None)
CLIP_BBOX_LV95 = tuple(_clip) if _clip else None   # None = full extent

CODE_FIELD_CANDIDATES = ("BLW_Code", "blw_code", "lnf_code", "LNF_CODE", "code")

GROUPS = {
    "mown": list(CFG.landuse.mown),
    "grazed": list(CFG.landuse.grazed),
    "intensive": list(CFG.landuse.intensive),
    "extensive": list(CFG.landuse.extensive),
}

# Present only in the latest delivery: a schema change, not land-use change,
# so exempt from persistence rather than failing it.
GRANDFATHER_CODES = [930]
