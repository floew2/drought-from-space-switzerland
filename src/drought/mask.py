"""Grassland mask — polygon-first.

WHAT CHANGED AND WHY. The previous version applied persistence per CELL, by
multiplying the purity raster with a binary "was grassland every year" layer.
Two things went wrong with that:

  * the analysed sample was a set of surviving CELLS belonging to no defined
    set of fields, so the number of parcels behind any figure was unknowable
    and the area accounting never resolved;
  * the multiplication zeroed the PSF halo outside parcels, giving an observed
    area SMALLER than the polygon area at zero threshold — physically
    impossible, since the sensor response necessarily spreads beyond every
    boundary.

The design here follows from the question being asked:

  1. MASTER PARCEL SET — parcels that are grassland in every year AND hold the
     same management category in every year. A fixed, countable set of fields.
     Persistence ends here and is never applied to a raster.
  2. FOUR CATEGORIES — by LNF code within that master set.
  3. PSF PURITY per category, from those polygons, on the 10 m Sentinel grid.

Category persistence matters because a parcel that was extensive until 2025 and
intensive now would otherwise carry an extensive management history into an
intensive baseline. The share excluded is reported, so it can be stated.

Output: four float32 bands of purity 0-1. No persistence bands.
"""

from __future__ import annotations

import re
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw

from .config import CFG, mask_path, parcels_dir

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
SHP_DIR = None                       # None -> data/nutzungen, searched
                                     # RECURSIVELY (one folder per year is fine)
PATTERN = "nutzungen_*.shp"
CRS = "EPSG:2056"                    # CH1903+ / LV95

RES = CFG.mask.res
# Polygons below this are cadastral offcuts, dropped so they cannot corrupt the
# overlap matching. It is NOT an observability filter: purity decides that.
SLIVER_FLOOR_M2 = CFG.mask.sliver_floor_m2
SUBPIXEL = CFG.mask.subpixel
PSF_MTF = CFG.mask.psf_mtf_nyquist
PSF_MODEL = getattr(CFG.mask, "psf_model", "box_gauss")
PURITY_MIN = CFG.mask.purity_min
CATEGORY_PERSISTENCE = getattr(CFG.mask, "category_persistence", True)
# A current parcel counts as holding a category in an earlier year if at least
# this share of its area falls inside that year's dissolved grassland of that
# category. Overlap, not identity — outlines are redrawn, split and merged
# between deliveries, and no identity match survives that.
MIN_OVERLAP = getattr(CFG.mask, "min_overlap", 0.8)
# How many of the available years a parcel must have held its category in.
# 0 = all of them.
MIN_YEARS_SAME = getattr(CFG.mask, "min_years_same", 0)
# Simplify geometries before overlaying. MEASURED AND TURNED OFF: on 60,000
# 300-vertex polygons it took 38 s without and 69 s with a 1 m tolerance —
# simplify costs more than the intersection it saves. Left as a knob in case a
# future dataset behaves differently, but 0 is the right default.
SIMPLIFY_M = getattr(CFG.mask, "simplify_m", 0.0)

# Restrict the output grid. None = full extent of the parcel data.
_clip = getattr(CFG.mask, "clip_bbox_lv95", None)
CLIP_BBOX_LV95 = tuple(_clip) if _clip else None

CODE_FIELD_CANDIDATES = ("BLW_Code", "blw_code", "lnf_code", "LNF_CODE", "code")

LU = CFG.landuse
GROUPS = {
    "mown": list(LU.mown),
    "grazed": list(LU.grazed),
    "intensive": list(LU.intensive),
    "extensive": list(LU.extensive),
}

# Codes appearing only in the most recent year(s) are a schema change, not
# land-use change, and would fail persistence for the wrong reason.
GRANDFATHER_CODES = [930]


def psf_params(res: float | None = None, mtf: float | None = None,
               model: str | None = None) -> dict:
    """The modelled sensor response, derived from the system MTF at Nyquist.

    For a Gaussian, MTF(f) = exp(-2 pi^2 sigma^2 f^2) with f_Nyquist = 1/(2p).

    box_gauss  The detector footprint is a box one pixel wide, with its own MTF
               sinc(1/2) = 0.637 at Nyquist. Only the remainder of the system MTF
               belongs to the optics: 0.25 / 0.637 = 0.393, so sigma_optics =
               4.35 m. The combined response (box convolved with that Gaussian)
               has FWHM 12.7 m; the Gaussian alone would be 10.2 m.
    gauss      One Gaussian fitted to the whole system MTF: sigma 5.30 m, FWHM
               12.5 m. The earlier implementation, which built the mask behind
               the published area numbers; kept for comparison. As sampled
               (a 2 x 2 sub-cell mean, not a point) its effective response is
               ~12.95 m — see `measured_response_fwhm`.

    Args:
        res: Pixel size in metres (default: config).
        mtf: System MTF at Nyquist (default: config).
        model: "box_gauss" or "gauss" (default: config).

    Returns:
        dict with model, sigma_m (the Gaussian term), fwhm_m (of the whole
        response), mtf_system, mtf_detector, mtf_optics.

    Raises:
        ValueError: For an unknown model.
    """
    from scipy.optimize import brentq
    from scipy.special import erf

    p_ = RES if res is None else res
    m_ = PSF_MTF if mtf is None else mtf
    model = PSF_MODEL if model is None else model
    f_n = 1.0 / (2.0 * p_)
    if model == "box_gauss":
        mtf_det = float(np.sinc(p_ * f_n))           # numpy sinc is sin(pi x)/(pi x)
    elif model == "gauss":
        mtf_det = 1.0                                # detector folded into the fit
    else:
        raise ValueError(f"unknown psf_model {model!r}; use 'box_gauss' or 'gauss'")
    mtf_opt = m_ / mtf_det
    sigma = float(np.sqrt(-np.log(mtf_opt) / (2.0 * np.pi**2 * f_n**2)))
    if model == "box_gauss":
        def prof(x):
            s2 = sigma * np.sqrt(2.0)
            return 0.5 * (erf((x + p_ / 2) / s2) - erf((x - p_ / 2) / s2))
        half = brentq(lambda x: prof(x) / prof(0.0) - 0.5, 0.0, 10 * p_)
        fwhm = 2.0 * half
    else:
        fwhm = 2.0 * np.sqrt(2.0 * np.log(2.0)) * sigma
    return {"model": model, "sigma_m": sigma, "fwhm_m": float(fwhm),
            "mtf_system": float(m_), "mtf_detector": mtf_det,
            "mtf_optics": float(mtf_opt)}


def psf_sigma_m(res: float | None = None, mtf: float | None = None,
                model: str | None = None) -> float:
    """Sigma of the Gaussian term of the modelled response, in metres.

    Under the default box_gauss model this is the OPTICAL sigma (4.35 m), not
    the width of the whole response — see `psf_params` for that.
    """
    return psf_params(res, mtf, model)["sigma_m"]


def measured_response_fwhm(model: str | None = None, sigma_m: float | None = None) -> float:
    """FWHM of the response that `purity_psf` actually applies, measured.

    `psf_params` states what the model should be; this measures what the code
    does, by passing a line source (one 2.5 m sub-cell wide) through
    `purity_psf` at four sub-cell offsets and reading the profile at 2.5 m
    spacing. The line source itself adds its 2.5 m width, so the expected
    values are box(2.5) (x) response:

        box_gauss   box(10) (x) G(4.35)            12.67 m -> ~12.8 measured
        gauss       2x2 sub-cell mean (x) G(5.30)  ~12.95 m -> ~13.1 measured
                    (the four-point "centre" estimate spans 5 m, so the old
                    implementation was wider than its nominal 12.48 m)
        gauss, sigma 4.35 — swapping sigma without averaging the cell:
                                                   ~10.8 m -> ~11.0 measured

    Args:
        model: "box_gauss" or "gauss" (default: config).
        sigma_m: Override for the Gaussian sigma.

    Returns:
        Measured FWHM in metres.
    """
    from shapely.geometry import box as _box

    n, res, xs, vs = 21, float(RES), [], []
    for k in range(SUBPIXEL):
        x0 = 10 * res + k * res / SUBPIXEL
        p = purity_psf([_box(x0, 0, x0 + res / SUBPIXEL, n * res)],
                       (0.0, n * res, res), (n, n), model=model, sigma_m=sigma_m)
        for c in range(n):
            xs.append((c + 0.5) * res - (x0 + res / SUBPIXEL / 2))
            vs.append(float(p[n // 2, c]))
    o = np.argsort(xs)
    x, v = np.asarray(xs)[o], np.asarray(vs)[o]
    v = v / v.max()
    i0, i1 = np.where(v >= 0.5)[0][[0, -1]]
    left = np.interp(0.5, [v[i0 - 1], v[i0]], [x[i0 - 1], x[i0]])
    right = np.interp(0.5, [v[i1 + 1], v[i1]], [x[i1 + 1], x[i1]])
    return float(right - left)


def psf_caption(model: str | None = None) -> str:
    """One-line description of the response model, for figure captions."""
    p = psf_params(model=model)
    if p["model"] == "box_gauss":
        return (f"{RES:g} m detector box \u2297 Gaussian optics (\u03c3 {p['sigma_m']:.2f} m), "
                f"FWHM {p['fwhm_m']:.1f} m, from system MTF {p['mtf_system']:g} at Nyquist")
    return (f"Gaussian \u03c3 {p['sigma_m']:.2f} m, FWHM {p['fwhm_m']:.1f} m, "
            f"from system MTF {p['mtf_system']:g} at Nyquist")


# --------------------------------------------------------------------------- #
# Reading
# --------------------------------------------------------------------------- #
def _year_files(shp_dir=None) -> dict:
    d = Path(shp_dir or SHP_DIR) if (shp_dir or SHP_DIR) else parcels_dir()
    out = {}
    for f in sorted(d.rglob(PATTERN)):
        m = re.search(r"(20\d\d)", f.stem)
        if m:
            out.setdefault(int(m.group(1)), f)
    if not out:
        raise SystemExit(
            f"No files matching {PATTERN} under {d.resolve()} (searched "
            "recursively).\nPut the yearly Nutzungsflaechen shapefiles there — "
            "one folder per year is fine — or set mask.SHP_DIR.")
    return out


def load_year(path: Path) -> gpd.GeoDataFrame:
    g = gpd.read_file(path)
    if g.crs is None:
        g = g.set_crs(CRS)
    elif g.crs.to_string() != CRS:
        g = g.to_crs(CRS)
    field = next((c for c in CODE_FIELD_CANDIDATES if c in g.columns), None)
    if field is None:
        raise SystemExit(f"no land-use code field in {path.name}; "
                         f"columns: {list(g.columns)[:12]}")
    g = g.rename(columns={field: "_code"})
    g["_code"] = pd.to_numeric(g["_code"], errors="coerce")
    g = g.dropna(subset=["_code"])
    g["_code"] = g["_code"].astype(int)
    g = _repair(g)
    return g[g.geometry.notna() & (g.geometry.area > 0)]


def _repair(g: gpd.GeoDataFrame, label: str = "") -> gpd.GeoDataFrame:
    """Make every geometry valid, or drop it.

    buffer(0) fixes most self-intersections but not all: a ring that touches
    itself at a point can survive it and then throw a side-location conflict
    during intersection, several minutes into a run. make_valid handles those,
    and anything still broken afterwards is discarded with a count — one
    unusable parcel is not worth losing the overlay for.
    """
    from shapely.validation import make_valid

    bad = ~g.geometry.is_valid
    n = int(bad.sum())
    if n:
        fixed = g.loc[bad, "geometry"].buffer(0)
        still = ~fixed.is_valid
        if still.any():
            fixed.loc[still] = fixed.loc[still].apply(make_valid)
        g = g.copy()
        g.loc[bad, "geometry"] = fixed

    # make_valid can return lines or collections; keep only the polygon parts
    types = g.geometry.geom_type
    odd = ~types.isin(["Polygon", "MultiPolygon"])
    if odd.any():
        from shapely.geometry import MultiPolygon
        def polys_only(geom):
            if geom is None or geom.is_empty:
                return None
            if geom.geom_type in ("Polygon", "MultiPolygon"):
                return geom
            parts = [p for p in getattr(geom, "geoms", [])
                     if p.geom_type in ("Polygon", "MultiPolygon")]
            return MultiPolygon(
                [q for p in parts
                 for q in (p.geoms if p.geom_type == "MultiPolygon" else [p])]
            ) if parts else None
        g = g.copy()
        g.loc[odd, "geometry"] = g.loc[odd, "geometry"].apply(polys_only)

    g = g[g.geometry.notna()]
    dropped = int((~g.geometry.is_valid).sum())
    if dropped:
        g = g[g.geometry.is_valid]
    if n:
        print(f"    repaired {n - dropped} invalid geometries"
              + (f", dropped {dropped} beyond repair" if dropped else "")
              + (f" [{label}]" if label else ""))
    return g


def category_of(code) -> str | None:
    """Management category for a code. None = not grassland.

    intensive and extensive are checked first because they are the analytical
    split; mown/grazed is the broader "is there a cut signal" split and a code
    can belong to both.
    """
    try:
        c = int(code)
    except (TypeError, ValueError):
        return None
    if c in GROUPS["intensive"]:
        return "intensive"
    if c in GROUPS["extensive"]:
        return "extensive"
    if c in GROUPS["grazed"]:
        return "grazed"
    if c in GROUPS["mown"]:
        return "mown_other"
    return None


# --------------------------------------------------------------------------- #
# 1. Master parcel set
# --------------------------------------------------------------------------- #
def _overlap_years(cur: gpd.GeoDataFrame, past: gpd.GeoDataFrame, codes,
                   min_overlap: float, simplify_m: float = None) -> pd.Series:
    """For each current parcel, what share of it does `past` cover with these codes?

    A SPATIAL JOIN, not a dissolve: union_all() on ~100,000 polygons is one
    enormous geometry that then has to be intersected against every parcel in
    Python. sjoin uses the index and returns only the pairs that touch.

    Two further shortcuts, because the pairwise intersection of real cadastral
    polygons is expensive and there are ~110,000 of them:

      CONTAINED  a parcel wholly inside a past polygon is covered 100% by
                 definition — no intersection needed. This is most of them.
    Simplifying first was tried and made it slower — see SIMPLIFY_M.

    Summing pairwise intersections is exact as long as the past polygons do not
    overlap one another, which parcel data does not.
    """
    sub = past[past["_code"].isin(codes)]
    if sub.empty:
        return pd.Series(0.0, index=cur.index)

    tol = SIMPLIFY_M if simplify_m is None else simplify_m
    left_geom = cur.geometry.simplify(tol) if tol else cur.geometry
    right_geom = sub.geometry.simplify(tol) if tol else sub.geometry

    # Build from the frame itself, not from a dict. Passing {"_i": cur.index}
    # with a separately-indexed GeoSeries makes pandas construct a fresh
    # RangeIndex for the column while the geometry keeps its own — so whenever
    # `cur` is a SUBSET (any per-category call) the two do not align and rows
    # are silently lost. That produced ~36% matches where 100% was correct.
    left = cur[["geometry"]].copy()
    left["geometry"] = left_geom
    left["_i"] = left.index
    right = sub[["geometry"]].copy()
    right["geometry"] = right_geom

    pairs = gpd.sjoin(left, right, how="inner", predicate="intersects")
    if pairs.empty:
        return pd.Series(0.0, index=cur.index)

    rg = gpd.GeoSeries(right.geometry.loc[pairs["index_right"]].to_numpy(),
                       index=pairs.index, crs=sub.crs)

    # Contained parcels contribute their whole area without any intersection.
    #
    # Every geometric predicate is guarded, because GEOS raises on malformed
    # rings that survive both buffer(0) and make_valid. But the guard is applied
    # PER CHUNK, not per pair: an earlier version fell back to a Python loop
    # over all ~250,000 pairs the moment one of them failed, which turned a
    # thirty-second overlay into an hour. Now one bad geometry costs one chunk.
    CH = 5000
    area = pd.Series(0.0, index=pairs.index)
    n_bad = 0
    for start in range(0, len(pairs), CH):
        sl = slice(start, min(start + CH, len(pairs)))
        gl = pairs.geometry.iloc[sl]
        gr = rg.iloc[sl]
        try:
            ins = gl.within(gr)
            a = pd.Series(0.0, index=gl.index)
            a[ins] = gl[ins].area
            rest = ~ins
            if rest.any():
                a[rest] = gl[rest].intersection(gr[rest]).area
        except Exception:
            a = pd.Series(0.0, index=gl.index)
            for ix, (x1, x2) in enumerate(zip(gl, gr)):
                try:
                    a.iloc[ix] = (x1.area if x1.within(x2)
                                  else x1.intersection(x2).area)
                except Exception:
                    try:
                        a.iloc[ix] = x1.buffer(0.01).intersection(
                            x2.buffer(0.01)).area
                    except Exception:
                        n_bad += 1
        area.iloc[sl] = a.to_numpy()
    if n_bad:
        print(f"    [!] {n_bad} pairs unusable even after nudging")

    got = area.groupby(pairs["_i"]).sum()
    frac = (got / cur.geometry.area.reindex(got.index)).reindex(cur.index)
    return frac.fillna(0.0).clip(upper=1.0)


def master_parcels(shp_dir=None, verbose=True, min_overlap=None,
                   min_years=None, sliver_floor_m2=None) -> gpd.GeoDataFrame:
    """Parcels that were grassland, and the same category, in earlier years.

    MATCHED BY OVERLAP, NOT BY IDENTITY. Parcel outlines are redrawn between
    yearly deliveries — boundaries adjusted, and a field split into two or two
    merged into one. Any method that tries to follow an individual parcel
    through the years fails on that: a split moves both centroids far from the
    original, and no matching tolerance recovers it.

    So we never identify parcels across years at all. For each year the
    grassland polygons are dissolved into one geometry per category, and a
    current parcel counts as "held that category in year Y" if at least
    `min_overlap` of its area falls inside that year's dissolved geometry.
    Splits, merges and redrawn edges all pass, because the question is about
    ground cover, not about parcel identity.

    A parcel qualifies for the master set if it held its current category in at
    least `min_years` of the available years. Codes in GRANDFATHER_CODES are
    exempt — they appear only in the most recent delivery, so there is nothing
    to check them against.

    THERE IS NO SIZE FILTER. Polygons below `sliver_floor_m2` (0.05 ha) are
    dropped as cadastral offcuts, in every year, so they cannot corrupt the
    overlap matching. Whether a field is large enough for Sentinel-2 is left to
    pixel purity, which already removes a field that yields no clean cell. An
    earlier 0.25 ha filter duplicated that job and, applied to past deliveries,
    failed stable fields whose redrawn outline dipped below it in one year.

    Args:
        shp_dir: Folder with the yearly shapefiles (default: config).
        verbose: Print the funnel.
        min_overlap: Share of a field that must lie inside a year's category.
        min_years: Years the category must be held (0/None = all).
        sliver_floor_m2: Override the sliver floor; 2500 reproduces the old
            0.25 ha behaviour.

    Returns:
        The master parcel set, EPSG:2056, with `_code` and `_cat`.
    """

    import time as _time
    files = _year_files(shp_dir)
    years = sorted(files)
    cur_year = years[-1]
    min_overlap = MIN_OVERLAP if min_overlap is None else min_overlap
    min_years = (MIN_YEARS_SAME if min_years is None else min_years) or len(years)

    floor = SLIVER_FLOOR_M2 if sliver_floor_m2 is None else sliver_floor_m2
    cur = load_year(files[cur_year])
    cur = cur[cur.geometry.area >= floor].copy()
    cur["_cat"] = cur["_code"].map(category_of)
    cur = cur[cur["_cat"].notna()].copy()
    n0 = len(cur)
    if verbose:
        print(f"  {n0:,} grassland parcels in {cur_year} "
              f"(slivers below {floor / 1e4:g} ha dropped)")

    cats = sorted(cur["_cat"].unique())
    codes_of = {"intensive": GROUPS["intensive"], "extensive": GROUPS["extensive"],
                "grazed": GROUPS["grazed"], "mown_other": GROUPS["mown"]}

    n_same = pd.Series(0, index=cur.index)
    n_grass = pd.Series(0, index=cur.index)
    grass_codes = sorted(set(GROUPS["mown"]) | set(GROUPS["grazed"]))

    for y in years[:-1]:
        t0 = _time.time()
        if verbose:
            print(f"  overlaying {y}...", end="", flush=True)
        g = load_year(files[y])
        g = g[g.geometry.area >= floor]

        n_grass += (_overlap_years(cur, g, grass_codes, min_overlap)
                    >= min_overlap).astype(int)
        for c in cats:
            sel = cur["_cat"] == c
            if not sel.any():
                continue
            frac = _overlap_years(cur[sel], g, codes_of[c], min_overlap)
            n_same.loc[sel] += (frac >= min_overlap).astype(int)
        if verbose:
            print(f" {_time.time() - t0:.0f} s")

    if cur_year in years:
        n_same += 1          # the current year holds its own category
        n_grass += 1

    grandfathered = cur["_code"].isin(GRANDFATHER_CODES)
    keep = (n_same >= min_years) | grandfathered

    out = cur[keep].drop(columns=["_cat"], errors="ignore").copy()
    out["_cat"] = cur.loc[keep, "_cat"]
    if verbose:
        lost_grass = int(((n_grass < min_years) & ~grandfathered).sum())
        lost_cat = int(((n_grass >= min_years) & (n_same < min_years)
                        & ~grandfathered).sum())
        print(f"\n  MASTER PARCEL SET ({cur_year})")
        print(f"    rule: held its category in >= {min_years} of {len(years)} "
              f"years,\n          matched by >= {min_overlap:.0%} area overlap")
        print(f"    {lost_grass:>8,}  dropped — not grassland in enough years")
        print(f"    {lost_cat:>8,}  dropped — grassland, but a different "
              "category")
        print(f"    {int(grandfathered.sum()):>8,}  exempt "
              f"(codes {GRANDFATHER_CODES}, present only in the latest year)")
        print(f"    {len(out):>8,}  remain "
              f"({100 * len(out) / max(n0, 1):.1f}%), "
              f"{out.geometry.area.sum() / 1e4:,.0f} ha")
        print("\n    per category:")
        for name, codes in GROUPS.items():
            sub = out[out["_code"].isin(codes)]
            print(f"      {name:<12}{len(sub):>8,} parcels"
                  f"{sub.geometry.area.sum() / 1e4:>12,.0f} ha")
        print("\n    years the category was held, distribution:")
        vc = n_same[~grandfathered].value_counts().sort_index()
        for k, v in vc.items():
            print(f"      {k:>2} of {len(years)}{v:>10,} parcels"
                  f"{100 * v / max(len(cur), 1):>8.1f}%")
    return out



# --------------------------------------------------------------------------- #
# 2. Rasterising and PSF purity
# --------------------------------------------------------------------------- #
def _burn(geoms, grid, shape) -> np.ndarray:
    """Rasterise polygons: 1 where a cell centre falls inside a polygon.

    Uses GDAL's rasteriser via rasterio, whose centre-in-polygon rule does not
    depend on where the window starts — so a tiled build gives exactly the
    single-pass result. PIL's polygon fill, used earlier to sidestep a
    rasterio/affine version conflict, flips edge pixels depending on the window
    offset (about 50 edge pixels in a 172 x 172 window), which made tiled purity
    differ from untiled purity by up to 0.15. It remains only as a fallback.

    Args:
        geoms: Iterable of shapely (Multi)Polygons.
        grid: (minx, maxy, res) of the raster.
        shape: (height, width) in cells.

    Returns:
        uint8 array of 0/1.
    """
    geoms = [g for g in geoms if g is not None and not g.is_empty]
    if not geoms:
        return np.zeros(shape, np.uint8)
    minx, maxy, res = grid
    try:
        from rasterio import features
        from rasterio.transform import from_origin
    except ImportError:                                   # pragma: no cover
        return _burn_pil(geoms, grid, shape)
    return features.rasterize(((g, 1) for g in geoms), out_shape=shape,
                              transform=from_origin(minx, maxy, res, res),
                              fill=0, dtype="uint8")


def _burn_pil(geoms, grid, shape) -> np.ndarray:
    """Fallback rasteriser (PIL). Edge pixels depend on the window offset."""
    minx, maxy, res = grid
    img = Image.new("1", (shape[1], shape[0]), 0)
    draw = ImageDraw.Draw(img)
    for geom in geoms:
        if geom is None or geom.is_empty:
            continue
        parts = geom.geoms if geom.geom_type == "MultiPolygon" else [geom]
        for poly in parts:
            xs, ys = poly.exterior.coords.xy
            px = [(x - minx) / res for x in xs]
            py = [(maxy - y) / res for y in ys]
            if len(px) > 2:
                draw.polygon(list(zip(px, py)), fill=1)
            for ring in poly.interiors:
                xs, ys = ring.coords.xy
                px = [(x - minx) / res for x in xs]
                py = [(maxy - y) / res for y in ys]
                if len(px) > 2:
                    draw.polygon(list(zip(px, py)), fill=0)
    return np.array(img, dtype=np.uint8)


def purity_psf(geoms, grid, shape, sub_n=None, sigma_m=None, model=None,
               tile: int = 512) -> np.ndarray:
    """PSF-weighted purity: the share of a cell's RESPONSE from inside a parcel.

    Rasterise at res/sub_n (2.5 m), convolve with the modelled response, and
    read the value for each 10 m cell. What contaminates the reflectance is the
    response, not the area, so this is the quantity a plain area fraction only
    approximates.

    Under box_gauss the Gaussian (optics) is applied on the fine grid and the
    detector box is applied by averaging the sub_n x sub_n sub-cells of each
    cell — which is exactly the box convolution evaluated at the cell centre.
    Under gauss the single Gaussian is sampled at the cell centre.

    Computed in tiles with a halo wider than the response, so memory stays
    bounded for a canton-sized grid and the result is identical to one pass.

    Args:
        geoms: Iterable of shapely (Multi)Polygons in the grid CRS.
        grid: (minx, maxy, res) of the output grid.
        shape: (height, width) of the output grid in cells.
        sub_n: Sub-cells per side (default: config, 4 -> 2.5 m).
        sigma_m: Override for the Gaussian sigma in metres.
        model: "box_gauss" or "gauss" (default: config).
        tile: Output tile size in cells.

    Returns:
        float32 array of purity in [0, 1], shape `shape`.
    """
    from scipy.ndimage import gaussian_filter
    from shapely.geometry import box as _box
    from shapely.strtree import STRtree

    sub_n = SUBPIXEL if sub_n is None else sub_n
    model = PSF_MODEL if model is None else model
    minx, maxy, res = grid
    height, width = shape
    fine_res = res / sub_n
    sigma = psf_sigma_m(model=model) if sigma_m is None else sigma_m
    sigma_px = sigma / fine_res
    # gaussian_filter truncates at 4 sigma; the halo must cover that plus the
    # half-cell of the detector box, so tile edges never see a missing parcel.
    halo = int(np.ceil((4.0 * sigma_px + 1) / sub_n)) + 1

    geoms = [g for g in geoms if g is not None and not g.is_empty]
    out = np.zeros((height, width), np.float32)
    if not geoms:
        return out
    tree = STRtree(geoms)

    for r0 in range(0, height, tile):
        r1 = min(r0 + tile, height)
        for c0 in range(0, width, tile):
            c1 = min(c0 + tile, width)
            R0, R1, C0, C1 = r0 - halo, r1 + halo, c0 - halo, c1 + halo
            x0, x1 = minx + C0 * res, minx + C1 * res
            y1, y0 = maxy - R0 * res, maxy - R1 * res
            hit = tree.query(_box(x0, y0, x1, y1))
            if not len(hit):
                continue
            nr, nc = R1 - R0, C1 - C0
            a = _burn([geoms[i] for i in hit], (x0, y1, fine_res),
                      (nr * sub_n, nc * sub_n)).astype(np.float32)
            a = gaussian_filter(a, sigma_px, mode="constant", cval=0.0)
            if model == "box_gauss":
                cell = a.reshape(nr, sub_n, nc, sub_n).mean(axis=(1, 3))
            elif sub_n % 2 == 0:             # the centre falls between sub-cells
                c = sub_n // 2
                cell = 0.25 * (a[c - 1::sub_n, c - 1::sub_n] + a[c - 1::sub_n, c::sub_n]
                               + a[c::sub_n, c - 1::sub_n] + a[c::sub_n, c::sub_n])
            else:
                cell = a[sub_n // 2::sub_n, sub_n // 2::sub_n]
            out[r0:r1, c0:c1] = cell[halo:halo + (r1 - r0), halo:halo + (c1 - c0)]
    # float32 convolution leaves a fully interior cell at 0.99999994, which
    # would make "purity >= 1" select nothing. Six decimals is far below any
    # threshold anyone sets.
    return np.clip(np.round(out, 6), 0.0, 1.0).astype(np.float32)


# --------------------------------------------------------------------------- #
# 3. Build
# --------------------------------------------------------------------------- #
def output_grid(master, verbose: bool = True) -> tuple[tuple, tuple]:
    """The 10 m output grid for a parcel set: ((minx, maxy, res), (height, width))."""
    b = master.total_bounds
    if CLIP_BBOX_LV95:
        x0 = max(b[0], CLIP_BBOX_LV95[0]); y0 = max(b[1], CLIP_BBOX_LV95[1])
        x1 = min(b[2], CLIP_BBOX_LV95[2]); y1 = min(b[3], CLIP_BBOX_LV95[3])
        if verbose:
            print("\n  grid clipped to mask.clip_bbox_lv95 in config.yaml "
                  "(null for the full extent)")
    else:
        x0, y0, x1, y1 = b
    pad = 2 * psf_params()["fwhm_m"]   # so the halo is not cut at the edge
    minx = np.floor((x0 - pad) / RES) * RES
    maxy = np.ceil((y1 + pad) / RES) * RES
    width = int((np.ceil((x1 + pad) / RES) * RES - minx) / RES)
    height = int((maxy - np.floor((y0 - pad) / RES) * RES) / RES)
    return (minx, maxy, RES), (height, width)


def psf_model_check(master, thresholds=(0.0, 0.5, 0.9), verbose: bool = True) -> pd.DataFrame:
    """Does the response model change the answer? Area ratios under both models.

    Builds purity for every category under box_gauss (the settled model) and
    gauss (the earlier single fitted Gaussian), on the same parcels and grid,
    and reports observed / polygon area at each threshold. The two responses
    have nearly the same second moment (5.22 m against 5.30 m), so the
    expectation is a small difference — this measures it instead of assuming it.

    Args:
        master: The master parcel set (from `master_parcels`).
        thresholds: Purity thresholds; 0 means "> 0".
        verbose: Print the table.

    Returns:
        DataFrame indexed by (category, threshold) with one ratio column per
        model and their difference.
    """
    grid, shape = output_grid(master, verbose=False)
    cell_ha = RES * RES / 1e4
    rows = []
    for name, codes in GROUPS.items():
        sub = master[master["_code"].isin(codes)]
        if sub.empty:
            continue
        poly_ha = float(sub.geometry.area.sum() / 1e4)
        bands = {m: purity_psf(sub.geometry, grid, shape, model=m)
                 for m in ("box_gauss", "gauss")}
        for t in thresholds:
            r = {m: float((p >= max(t, 1e-9)).sum()) * cell_ha / poly_ha
                 for m, p in bands.items()}
            rows.append({"category": name, "threshold": t, **r,
                         "difference": r["box_gauss"] - r["gauss"]})
    out = pd.DataFrame(rows).set_index(["category", "threshold"])
    if verbose:
        print("  observed / polygon area, by response model:")
        print(out.round(3).to_string())
    return out


def main(shp_dir=None, out_tif=None):
    master = master_parcels(shp_dir)
    grid, (height, width) = output_grid(master)

    print(f"  grid {width} x {height} px @ {RES} m, {CRS}")
    print(f"  response: {psf_caption()}; fine grid {RES / SUBPIXEL:g} m")

    cell_ha = (RES * RES) / 1e4
    bands = []
    print(f"\n  {'category':<12}{'parcels':>9}{'polygon ha':>12}"
          f"{'ha at >0':>11}{'ratio':>7}{'ha at >=' + f'{PURITY_MIN:g}':>13}"
          f"{'ratio':>7}")
    for name, codes in GROUPS.items():
        sub = master[master["_code"].isin(codes)]
        if sub.empty:
            bands.append((np.zeros((height, width), np.float32), name))
            print(f"  {name:<12}{'—':>9}")
            continue
        p = purity_psf(sub.geometry, grid, (height, width))
        poly_ha = float(sub.geometry.area.sum() / 1e4)
        h0 = float((p > 0).sum()) * cell_ha
        h9 = float((p >= PURITY_MIN).sum()) * cell_ha
        print(f"  {name:<12}{len(sub):>9,}{poly_ha:>12,.0f}{h0:>11,.0f}"
              f"{h0 / poly_ha:>7.2f}{h9:>13,.0f}{h9 / poly_ha:>7.2f}")
        bands.append((p, name))

    print("\n  the ratio at >0 MUST exceed 1: the PSF response spreads beyond "
          "every\n  boundary. If it does not, something is wrong.")

    out_tif = Path(out_tif) if out_tif else mask_path()
    _write(bands, grid, (height, width), out_tif)
    print(f"\n  wrote {out_tif} ({out_tif.stat().st_size / 1e6:.1f} MB)")
    print("  band order: " + ", ".join(f"b{i} {n}"
                                       for i, (_, n) in enumerate(bands, 1)))
    print("\n  FLOAT purity 0-1. No persistence multiplication — persistence is "
          "already in\n  the parcel selection. Threshold in GEE with .gte(0.9), "
          "and upload with\n  pyramidingPolicy {'.default': 'sample'} or reduce "
          "at the native 10 m.")
    return master


def _write(bands, grid, shape, out_tif):
    minx, maxy, res = grid
    height, width = shape
    try:
        import rasterio
        from rasterio.transform import from_origin
        with rasterio.open(
            out_tif, "w", driver="GTiff", height=height, width=width,
            count=len(bands), dtype="float32", crs=CRS,
            transform=from_origin(minx, maxy, res, res),
            compress="deflate", predictor=3, tiled=True
        ) as ds:
            for i, (arr, name) in enumerate(bands, 1):
                ds.write(arr.astype(np.float32), i)
                ds.set_band_description(i, name)
        print("    written with rasterio")
    except Exception as exc:
        from osgeo import gdal, osr
        print(f"    rasterio unavailable ({exc}); using GDAL")
        drv = gdal.GetDriverByName("GTiff")
        ds = drv.Create(str(out_tif), width, height, len(bands),
                        gdal.GDT_Float32, ["COMPRESS=DEFLATE", "TILED=YES"])
        ds.SetGeoTransform((minx, res, 0, maxy, 0, -res))
        srs = osr.SpatialReference(); srs.ImportFromEPSG(2056)
        ds.SetProjection(srs.ExportToWkt())
        for i, (arr, name) in enumerate(bands, 1):
            band = ds.GetRasterBand(i)
            band.WriteArray(arr.astype(np.float32))
            band.SetDescription(name)
        ds = None
def _site_bounds_lv95(site):
    """Site bbox from config, transformed WGS84 -> LV95."""
    from pyproj import Transformer

    from .config import site_cfg
    b = site_cfg(site).bbox_wgs84
    t = Transformer.from_crs("EPSG:4326", CRS, always_xy=True)
    x0, y0 = t.transform(b[0], b[1])
    x1, y1 = t.transform(b[2], b[3])
    return min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1)


def figure_purity_map(site, group="mown", out=None, inset=None,
                     inset_frac=0.14, threshold=None):
    """Continuous purity for ONE landscape, with a zoomed inset.

    Deliberately NOT thresholded: the gradient is the point. A cell at 0.5 looks
    different from one at 0.95, and that difference is exactly what a binary
    mask throws away. The threshold appears only as a marker on the colour bar,
    so the reader can see where it sits within the distribution rather than
    being shown its result.

    `inset` is (x0, y0, x1, y1) in LV95. If omitted, a window of `inset_frac`
    of the extent is taken from the densest part of the mask, so the parcel-level
    texture is visible rather than a mostly empty corner.
    """
    import matplotlib.pyplot as plt
    import rasterio
    from matplotlib.patches import Rectangle

    from .config import figures_dir

    threshold = PURITY_MIN if threshold is None else threshold
    band = list(GROUPS).index(group) + 1
    x0, y0, x1, y1 = _site_bounds_lv95(site)
    with rasterio.open(mask_path()) as ds:
        win = rasterio.windows.from_bounds(x0, y0, x1, y1, ds.transform)
        a = ds.read(band, window=win)
        tr = ds.window_transform(win)
    if not a.size or not (a > 0).any():
        print(f"no {group} parcels in {site}")
        return None

    res = tr.a
    # Crop to where this group actually occurs, otherwise a landscape whose
    # parcels sit in one corner renders as mostly empty black.
    rows, cols = np.where(a > 0)
    pad = 8
    r0, r1 = max(rows.min() - pad, 0), min(rows.max() + pad + 1, a.shape[0])
    c0, c1 = max(cols.min() - pad, 0), min(cols.max() + pad + 1, a.shape[1])
    a = a[r0:r1, c0:c1]
    ox, oy = tr.c + c0 * res, tr.f - r0 * res
    ext = [ox, ox + a.shape[1] * res, oy - a.shape[0] * res, oy]

    # Densest window, so the inset shows fields rather than empty ground.
    if inset is None:
        n = max(int(min(a.shape) * inset_frac), 12)
        occ = (a > 0).astype(float)
        cs = occ.cumsum(0).cumsum(1)
        best, br, bc = -1, 0, 0
        for r in range(0, a.shape[0] - n, max(n // 3, 1)):
            for c in range(0, a.shape[1] - n, max(n // 3, 1)):
                tot = (cs[r + n, c + n] - cs[r, c + n] - cs[r + n, c] + cs[r, c]
                       if r and c else occ[r:r + n, c:c + n].sum())
                if tot > best:
                    best, br, bc = tot, r, c
        inset = (ext[0] + bc * res, ext[3] - (br + n) * res,
                 ext[0] + (bc + n) * res, ext[3] - br * res)
    ir0 = int((ext[3] - inset[3]) / res); ic0 = int((inset[0] - ext[0]) / res)
    ir1 = int((ext[3] - inset[1]) / res); ic1 = int((inset[2] - ext[0]) / res)
    sub = a[ir0:ir1, ic0:ic1]

    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(11, 8.4), dpi=180, constrained_layout=True)
    m = np.ma.masked_where(a <= 0, a)
    # Black -> white, so purity reads as brightness: pure cells glow, mixed
    # edges fade out. A hue ramp made the rim a colour rather than a gradient.
    cmap = plt.cm.gray.copy(); cmap.set_bad("black")
    im = ax.imshow(m, extent=ext, cmap=cmap, vmin=0, vmax=1,
                   interpolation="nearest")
    ax.add_patch(Rectangle((inset[0], inset[1]), inset[2] - inset[0],
                           inset[3] - inset[1], fill=False, edgecolor="#58a6ff",
                           lw=1.8, zorder=5))

    # --- inset, upper right
    axi = ax.inset_axes([0.615, 0.615, 0.375, 0.375])
    mi = np.ma.masked_where(sub <= 0, sub)
    axi.imshow(mi, cmap=cmap, vmin=0, vmax=1, interpolation="nearest")
    axi.set_xticks([]); axi.set_yticks([])
    for sp in axi.spines.values():
        sp.set_color("#58a6ff"); sp.set_linewidth(1.8)
    axi.text(0.02, 0.975, f"{(inset[2]-inset[0])/1000:.1f} km wide",
             transform=axi.transAxes, fontsize=8, color="#58a6ff", va="top",
             bbox=dict(facecolor="black", edgecolor="none", pad=2))

    # --- scale bar, lower left
    span = ext[1] - ext[0]
    bar = max(round(span / 5 / 1000), 1) * 1000
    bx, by = ext[0] + span * 0.04, ext[2] + (ext[3] - ext[2]) * 0.055
    ax.plot([bx, bx + bar], [by, by], color="white", lw=3, solid_capstyle="butt",
            zorder=6, path_effects=[])
    ax.text(bx + bar / 2, by + (ext[3] - ext[2]) * 0.018, f"{bar / 1000:g} km",
            ha="center", fontsize=9, color="white", zorder=6,
            bbox=dict(facecolor="black", edgecolor="none", pad=1.5))

    # --- north arrow, lower right
    nx, ny = ext[1] - span * 0.055, ext[2] + (ext[3] - ext[2]) * 0.055
    ax.annotate("", xy=(nx, ny + (ext[3] - ext[2]) * 0.075), xytext=(nx, ny),
                zorder=6,
                arrowprops=dict(arrowstyle="-|>", color="white", lw=1.6))
    ax.text(nx, ny - (ext[3] - ext[2]) * 0.012, "N", ha="center", va="top",
            fontsize=11, color="white", fontweight="bold", zorder=6,
            bbox=dict(facecolor="black", edgecolor="none", pad=1.5))

    # --- LV95 coordinates, sparse
    ax.set_xlabel("E [m], CH1903+ / LV95", fontsize=9, color="#8b949e")
    ax.set_ylabel("N [m]", fontsize=9, color="#8b949e")
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    ax.ticklabel_format(style="plain", useOffset=False)
    ax.tick_params(colors="#8b949e", labelsize=8)
    for sp in ax.spines.values():
        sp.set_color("#484f58")

    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.015)
    cb.set_label(f"purity — share of the {RES:g} m cell's response from inside a parcel",
                 fontsize=9, color="white")
    cb.ax.axhline(threshold, color="#ff3b30", lw=2)
    cb.ax.tick_params(colors="white", labelsize=8)
    cb.ax.text(1.6, threshold, f" {threshold:.2f}", color="#ff3b30",
               fontsize=8.5, va="center", transform=cb.ax.get_yaxis_transform())

    kept = 100 * (a[a > 0] >= threshold).mean()
    ax.set_title(f"Pixel purity  ·  {site}  ·  {group}", fontsize=12,
                 fontweight="bold", color="#ff3b30", loc="left")
    fig.supxlabel(f"Continuous, not thresholded — the gradient is the point. "
                  f"Black = no parcel of this group. At {threshold:.2f}, "
                  f"{kept:.0f}% of cells are retained.",
                  fontsize=7.5, color="#6e7681")
    out = out or (figures_dir("mask") / f"purity_map_{site}_{group}.png")
    fig.savefig(out, dpi=180, facecolor=fig.get_facecolor())
    print(f"  -> {out}")
    plt.show()
    return a

# --------------------------------------------------------------------------- #
# Diagnostics
# --------------------------------------------------------------------------- #
def area_check(site=None, thresholds=(0.0, 0.25, 0.5, 0.75, 0.9), out=None,
               master=None):
    """Observed area against polygon area, per category, across thresholds.

    Now a meaningful comparison, because both sides describe the SAME fixed set
    of parcels — the master set. Previously the raster carried a cell-level
    persistence filter the polygon reference did not, so the two measured
    different things and the ratio could fall below 1 at zero threshold, which
    is impossible.

    Read the ratio at >0 as a check: it must exceed 1, because the PSF response
    spreads beyond every boundary. Crossing 1 is NOT the target — that happens
    where roughly half the retained cells are contaminated. The threshold is
    justified by the anomaly being flat across the whole range.
    """
    import matplotlib.pyplot as plt
    import rasterio
    from shapely.geometry import box as _box

    from .config import figures_dir

    if master is None:
        master = master_parcels(verbose=False)

    # ALWAYS clip the polygon reference to the ground the raster covers. The
    # grid is limited by CLIP_BBOX_LV95, so parcels outside it exist in the
    # polygons but cannot exist in the raster — counting them makes the
    # reference too large and the ratio drops below 1 at zero threshold, which
    # is impossible. It bit the grazed category hardest, because Soemmerung
    # units are large and mostly outside the box.
    with rasterio.open(mask_path()) as ds:
        rb = ds.bounds
        if site:
            sx0, sy0, sx1, sy1 = _site_bounds_lv95(site)
            x0, y0 = max(rb.left, sx0), max(rb.bottom, sy0)
            x1, y1 = min(rb.right, sx1), min(rb.top, sy1)
        else:
            x0, y0, x1, y1 = rb.left, rb.bottom, rb.right, rb.top
        win = rasterio.windows.from_bounds(x0, y0, x1, y1, ds.transform)
        names = list(GROUPS)
        arrs = {n: ds.read(i, window=win) for i, n in enumerate(names, 1)}

    clip = _box(x0, y0, x1, y1)
    master = master[master.geometry.intersects(clip)].copy()
    master["geometry"] = master.geometry.intersection(clip)
    master = master[master.geometry.area > 0]
    cell_ha = (RES * RES) / 1e4

    out_d = {}
    print(f"{site or 'full extent'}  ·  master parcel set  ·  reference "
          f"clipped to the raster extent\n")
    print(f"  {'category':<12}{'parcels':>9}{'polygon ha':>12}{'median ha':>11}"
          f"{'shape':>8}")
    for n in names:
        sub = master[master["_code"].isin(GROUPS[n])]
        a_ha = sub.geometry.area / 1e4
        shp = (sub.geometry.length / np.sqrt(sub.geometry.area)).median() \
            if len(sub) else float("nan")
        obs = {t: float((arrs[n] >= max(t, 1e-9)).sum()) * cell_ha
               for t in thresholds}
        poly = float(a_ha.sum())
        out_d[n] = {"parcels": len(sub), "polygon_ha": poly,
                    "median_ha": float(a_ha.median()) if len(sub) else np.nan,
                    "median_shape": float(shp),
                    "observed_ha": obs,
                    "ratio": {t: (v / poly if poly else np.nan)
                              for t, v in obs.items()}}
        print(f"  {n:<12}{len(sub):>9,}{poly:>12,.0f}"
              f"{out_d[n]['median_ha']:>11.2f}{shp:>8.2f}")

    print("\n  observed area [ha] and ratio, by purity threshold:")
    print("  " + f"{'category':<12}" + "".join(f"{'>' + f'{t:g}':>16}"
                                               for t in thresholds))
    for n in names:
        cells = "".join(f"{out_d[n]['observed_ha'][t]:>9,.0f}"
                        f"{out_d[n]['ratio'][t]:>7.2f}" for t in thresholds)
        print(f"  {n:<12}{cells}")
    print("\n  the ratio at >0 must exceed 1 — the PSF spreads past every "
          "boundary.")

    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(8.0, 5.6), dpi=180, constrained_layout=True)
    col = {"intensive": "#3fb950", "extensive": "#58a6ff",
           "mown": "#d29922", "grazed": "#8250df"}
    grid_t = np.round(np.arange(0.0, 0.96, 0.05), 2)
    with rasterio.open(mask_path()) as ds:
        for i, n in enumerate(names, 1):
            a = ds.read(i, window=win)
            poly = out_d[n]["polygon_ha"]
            if not poly:
                continue
            r = [float((a >= max(t, 1e-9)).sum()) * cell_ha / poly
                 for t in grid_t]
            ax.plot(grid_t, r, "-o", ms=3.5, lw=2.2, color=col[n],
                    label=f"{n}  ({out_d[n]['median_ha']:.2f} ha, "
                          f"shape {out_d[n]['median_shape']:.1f})")
    ax.axhline(1.0, color="#8b949e", ls=":", lw=1.4)
    ax.axvline(PURITY_MIN, color="#ff3b30", ls="--", lw=1.8,
               label=f"threshold in use ({PURITY_MIN:.2f})")
    ax.set_xlabel("purity threshold", fontsize=10, color="white")
    ax.set_ylabel("observed area / polygon area", fontsize=10, color="white")
    ax.set_xlim(0, 0.97)
    ax.set_title(f"How much of the master parcel set do we measure?"
                 f"  ·  {site or 'full extent'}", fontsize=11.5,
                 fontweight="bold", color="#ff3b30", loc="left")
    ax.legend(fontsize=8.5, frameon=True, facecolor="#161b22",
              edgecolor="none", labelcolor="white", loc="upper right")
    ax.grid(alpha=0.22, ls=":", color="gray")
    ax.tick_params(colors="white", labelsize=9)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color("#484f58")
    fig.supxlabel("Above 1 at low thresholds is the PSF halo, and is required. "
                  "Crossing 1 is not the target — it happens where about half "
                  "the retained\ncells are contaminated. Legend gives median "
                  "parcel size and shape index (4.0 = a perfect square).",
                  fontsize=7, color="#6e7681")
    out = out or (figures_dir("mask") / f"area_check_{site or 'all'}.png")
    fig.savefig(out, dpi=180, facecolor=fig.get_facecolor())
    print(f"\n  -> {out}")
    plt.show()
    return out_d


# --------------------------------------------------------------------------- #
# Land-use change diagnostics
# --------------------------------------------------------------------------- #
def category_history(shp_dir=None, min_overlap=None, verbose=True,
                     sliver_floor_m2=None) -> pd.DataFrame:
    """How many years each current parcel held its category, by area overlap.

    The companion to master_parcels: same matching, but it reports the whole
    distribution instead of applying a cut-off, so the rule can be chosen from
    evidence.

    Reports separately how often a parcel was ANY grassland and how often it was
    ITS OWN category. The gap between the two is real management change; a low
    count on both means the ground simply was not grassland.
    """
    import time as _time
    files = _year_files(shp_dir)
    years = sorted(files)
    cur_year = years[-1]
    min_overlap = MIN_OVERLAP if min_overlap is None else min_overlap

    floor = SLIVER_FLOOR_M2 if sliver_floor_m2 is None else sliver_floor_m2
    cur = load_year(files[cur_year])
    cur = cur[cur.geometry.area >= floor].copy()
    cur["_cat"] = cur["_code"].map(category_of)
    cur = cur[cur["_cat"].notna()].copy()

    grass_codes = sorted(set(GROUPS["mown"]) | set(GROUPS["grazed"]))
    codes_of = {"intensive": GROUPS["intensive"], "extensive": GROUPS["extensive"],
                "grazed": GROUPS["grazed"], "mown_other": GROUPS["mown"]}
    cats = sorted(cur["_cat"].unique())

    hist = pd.DataFrame(index=cur.index)
    hist["area_ha"] = cur.geometry.area / 1e4
    hist["cat"] = cur["_cat"]
    n_grass = pd.Series(1, index=cur.index)      # the current year counts
    n_same = pd.Series(1, index=cur.index)

    for y in years[:-1]:
        t0 = _time.time()
        if verbose:
            print(f"  overlaying {y}...", end="", flush=True)
        g = load_year(files[y])
        g = g[g.geometry.area >= floor]
        n_grass += (_overlap_years(cur, g, grass_codes, min_overlap)
                    >= min_overlap).astype(int)
        for c in cats:
            sel = cur["_cat"] == c
            if not sel.any():
                continue
            n_same.loc[sel] += (
                _overlap_years(cur[sel], g, codes_of[c], min_overlap)
                >= min_overlap).astype(int)
        if verbose:
            print(f" {_time.time() - t0:.0f} s")

    hist["years_grassland"] = n_grass
    hist["years_same_cat"] = n_same
    hist["grandfathered"] = cur["_code"].isin(GRANDFATHER_CODES)

    if verbose:
        n = len(hist)
        ny = len(years)
        print(f"\n  {n:,} grassland parcels in {cur_year}, over {ny} years, "
              f"matched by >= {min_overlap:.0%} area overlap\n")
        print(f"  {'years':>7}{'was grassland':>16}{'held its category':>20}")
        for k in range(ny, 0, -1):
            a = int((hist.years_grassland == k).sum())
            b = int((hist.years_same_cat == k).sum())
            print(f"  {k:>5} y{a:>12,} {100*a/n:>4.0f}%{b:>13,} {100*b/n:>4.0f}%")
        print("\n  WHAT EACH RULE WOULD KEEP")
        for k in range(ny, max(ny - 5, 0), -1):
            r = ((hist.years_same_cat >= k) | hist.grandfathered).mean()
            print(f"    category held in >= {k} of {ny} years{100 * r:>14.1f}%")
        chg = int(((hist.years_grassland == ny) & (hist.years_same_cat < ny)).sum())
        print(f"\n  {chg:,} parcels ({100 * chg / n:.1f}%) were grassland in "
              "every year but\n  changed category at some point — that is real "
              "management change.")
    return hist


def matching_check(shp_dir=None, overlaps=(0.5, 0.7, 0.8, 0.9, 0.95)):
    """How sensitive is the master set to the overlap threshold?

    A gentle slope means the matching is robust — parcels are either clearly
    the same ground or clearly not. A steep one means many parcels sit near the
    boundary, and the threshold is doing more work than it should.
    """
    print(f"{'overlap':>9}{'held category in ALL years':>29}")
    out = {}
    for o in overlaps:
        h = category_history(shp_dir, min_overlap=o, verbose=False)
        ny = h.years_same_cat.max()
        r = 100 * ((h.years_same_cat >= ny) | h.grandfathered).mean()
        out[o] = r
        print(f"{o:>8.0%}{r:>28.1f}%")
    print("\n  a gentle slope means the matching is robust; a steep one means "
          "the\n  threshold is doing more work than it should.")
    return out


# --------------------------------------------------------------------------- #
# Explainer figure — why purity exists, what it looks like, what it costs
# --------------------------------------------------------------------------- #

def figure_purity_explainer(site="linthebene", group="intensive", out=None,
                            zoom_m=400, thresholds=None, master=None,
                            areas=None, zoom_levels=(0.0, 0.6, 1.0)):
    # zoom_levels sets the three panels in A, e.g. zoom_levels=(0.3, 0.6, 0.9).
    # Any count works; whichever equals PURITY_MIN is framed and labelled red.
    """Three panels that build the argument in order.

    The point-spread function that explains all of this has its own figure —
    see figure_psf().

    A  WHAT THAT LOOKS LIKE ON THE GROUND. The same patch of real mask at three
       thresholds. At 0 every cell the response touches is kept, including
       ground outside the parcels; at 1.0 only the cores survive, because under
       a point-spread function a cell reaches 1 only well inside a large field.
       Watching the fields erode inward makes the cost concrete in a way the
       curve in C cannot.

    B  WHAT IT COSTS. Measured area against polygon area as the threshold moves,
       for all four categories. The spread between them is parcel geometry:
       small, convoluted fields are mostly rim.
    """
    import matplotlib.pyplot as plt
    import rasterio
    from matplotlib.patches import Rectangle

    from .config import figures_dir

    thresholds = list(thresholds or np.round(np.arange(0.0, 0.96, 0.05), 2))
    names = list(GROUPS)

    # --- data -------------------------------------------------------------
    x0, y0, x1, y1 = _site_bounds_lv95(site)
    with rasterio.open(mask_path()) as ds:
        b_ = ds.bounds
        xa, ya = max(b_.left, x0), max(b_.bottom, y0)
        xb, yb = min(b_.right, x1), min(b_.top, y1)
        band = names.index(group) + 1
        win = rasterio.windows.from_bounds(xa, ya, xb, yb, ds.transform)
        big = ds.read(band, window=win)
        full = {n: ds.read(names.index(n) + 1,
                           window=rasterio.windows.from_bounds(
                               b_.left, b_.bottom, b_.right, b_.top,
                               ds.transform))
                for n in names}

    # densest window for the zoom
    npx = max(int(zoom_m / RES), 20)
    occ = (big > 0).astype(float)
    best, br, bc = -1, 0, 0
    step = max(npx // 3, 1)
    for r_ in range(0, max(big.shape[0] - npx, 1), step):
        for c_ in range(0, max(big.shape[1] - npx, 1), step):
            v = occ[r_:r_ + npx, c_:c_ + npx].sum()
            if v > best:
                best, br, bc = v, r_, c_
    zoom = big[br:br + npx, bc:bc + npx]

    # Panel C needs ONE number per category: the mapped polygon area. Nothing
    # else here touches polygons — panel A is the Gaussian, panel B reads the
    # raster.
    #
    # Those four numbers are printed by main(). Pass them and the figure never
    # goes near a shapefile:
    #     areas={"mown": 44628, "grazed": 36223,
    #            "intensive": 42034, "extensive": 2585}
    # Or pass master= to sum them from a set you already built. Only with
    # neither does it rebuild, and rebuilding means eight yearly overlays —
    # minutes to hours for four numbers you already have.
    cell_ha = (RES * RES) / 1e4
    if areas is None:
        if master is None:
            print("  building the parcel set just to total four areas — pass "
                  "areas= or master= to skip this")
            master = master_parcels(verbose=False)
        areas = {n: float(master[master["_code"].isin(GROUPS[n])]
                          .geometry.area.sum() / 1e4) for n in names}

    curves = {}
    for n in names:
        poly = float(areas.get(n, 0.0))
        if poly <= 0:
            continue
        curves[n] = [float((full[n] >= max(t, 1e-9)).sum()) * cell_ha / poly
                     for t in thresholds]

    # --- figure -----------------------------------------------------------
    plt.style.use("dark_background")
    nz = len(zoom_levels)
    # The zoom strip and the cost curve get equal width, so neither dominates.
    fig = plt.figure(figsize=(6.2 + 2.15 * nz, 5.0), dpi=180)
    gs = fig.add_gridspec(1, nz + 1,
                          width_ratios=[0.80] * nz + [1.45],
                          wspace=0.20, left=0.035, right=0.985,
                          top=0.76, bottom=0.205)

    # A — the same ground at several thresholds
    zoom_axes = []
    for zi, lev in enumerate(zoom_levels):
        ax = fig.add_subplot(gs[zi])
        zoom_axes.append(ax)
        keep = zoom >= max(lev, 1e-9)
        shown = np.where(keep, zoom, np.nan)
        ax.imshow(np.ma.masked_invalid(shown), cmap="gray", vmin=0, vmax=1,
                  interpolation="nearest")
        ax.set_xticks([]); ax.set_yticks([])
        pct = 100 * keep.sum() / max((zoom > 0).sum(), 1)
        col_ = "#ff3b30" if abs(lev - PURITY_MIN) < 1e-9 else "white"
        ax.set_title(f"\u2265 {lev:g}", fontsize=11, fontweight="bold",
                     color=col_, loc="left", pad=4)
        ax.set_xlabel(f"{pct:.0f}% of cells kept", fontsize=9.5, color=col_,
                      labelpad=3)
        for sp in ax.spines.values():
            sp.set_color("#ff3b30" if col_ == "#ff3b30" else "#484f58")
        if zi == 0:
            sb = int(100 / RES)
            ax.add_patch(Rectangle((1, npx - 7), sb + 2, 6, color="black",
                                   alpha=0.75))
            ax.add_patch(Rectangle((2, npx - 3.4), sb, 1.3, color="white"))
            ax.text(2 + sb / 2, npx - 4.2, "100 m", ha="center", va="bottom",
                    fontsize=7.5, color="white")
    fig.text(0.035, 0.845, "A   The same ground, three pixel purity thresholds",
             fontsize=12, fontweight="bold", color="white", ha="left")

    # C — what the threshold costs
    ax = fig.add_subplot(gs[-1])
    col = {"intensive": "#3fb950", "extensive": "#58a6ff",
           "mown": "#d29922", "grazed": "#8250df"}
    # mown and intensive differ by under half a percent — intensive is ~95% of
    # mown by area — so a solid line simply hides the other. The broader groups
    # are drawn dashed and wider, underneath, so both stay visible.
    style = {"mown": dict(ls="--", lw=3.0, marker="", alpha=0.75, zorder=1),
             "grazed": dict(ls="--", lw=3.0, marker="", alpha=0.75, zorder=1),
             "intensive": dict(ls="-", lw=2.0, marker="o", ms=3, zorder=3),
             "extensive": dict(ls="-", lw=2.0, marker="o", ms=3, zorder=3)}
    for n in ("mown", "grazed", "intensive", "extensive"):
        if n not in curves:
            continue
        ax.plot(thresholds, curves[n], color=col[n], label=n, **style[n])
    ax.axhline(1, color="#8b949e", ls=":", lw=1.3)
    ax.axvline(PURITY_MIN, color="#ff3b30", ls="--", lw=1.8)
    ax.annotate(f"{PURITY_MIN:g}", xy=(PURITY_MIN, 0), xycoords=("data",
                "axes fraction"), xytext=(0, -17), textcoords="offset points",
                ha="center", va="top", fontsize=9.5, fontweight="bold",
                color="#ff3b30")
    ax.set_xlim(0, 0.96)
    top = max(2.0, min(4.0, max(max(v) for v in curves.values()) * 1.08))
    ax.set_ylim(0, top)
    ax.set_xlabel("purity threshold", fontsize=9.5, color="white")
    ax.set_ylabel("measured area / mapped area", fontsize=9.5, color="white")
    ax.set_title("B   What the threshold costs", fontsize=12,
                 fontweight="bold", color="white", loc="left", pad=8)
    ax.legend(fontsize=8, frameon=True, facecolor="#161b22", edgecolor="none",
              labelcolor="white", loc="upper right", ncol=1,
              bbox_to_anchor=(0.985, 0.985), framealpha=0.92,
              borderpad=0.6)

    for a in (ax,):
        a.tick_params(colors="white", labelsize=8.5)
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            a.spines[sp].set_color("#484f58")
        a.grid(alpha=0.18, ls=":", color="gray")

    fig.suptitle("Mixed pixels: why the mask keeps only field interiors",
                 fontsize=13, fontweight="bold", color="#ff3b30",
                 x=0.035, ha="left", y=0.945)
    fig.text(0.035, 0.018,
             f"Response: {psf_caption()}. "
             f"Panel A: {site} at {RES} m, {zoom_m} m across. Above 1 in "
             "panel B is the response reaching past the parcel; below 1 is "
             "edges being dropped.\nDashed = the broader groups (mown, "
             "grazed). mown and intensive nearly coincide, because intensive "
             "is most of mown by area.",
             fontsize=7.5, color="#6e7681", ha="left", va="bottom")

    out = out or (figures_dir("mask") / f"purity_explainer_{site}_{group}.png")
    fig.savefig(out, dpi=180, facecolor=fig.get_facecolor())
    print(f"  -> {out}")
    plt.show()
    return curves


def figure_mask_overview(sites=None, thresholds=(0.5, 0.9), group="intensive",
                         zoom_m=500, rgb=None, out=None, master=None,
                         stretch=(0.02, 0.98), windows=None, labels=None,
                         rgb_bands=(1, 2, 3), verbose=True):
    """Four columns per site: imagery, polygons, and purity at two thresholds.

    Reads left to right as the mask is actually built — what the sensor saw,
    what the cadastre says, and what survives once the sensor's own blur is
    accounted for. Three rows let a reader see that the cost of the threshold
    is a property of the landscape, not of the method.

    Column 1 needs a Sentinel-2 GeoTIFF per site, exported from GEE in LV95 and
    placed at data/gee/<site>/rgb_<site>.tif (or passed in `rgb`). Three bands,
    conventionally B8/B4/B3 for false colour, where vegetation reads red. If
    the file is absent the panel says so rather than inventing something.

    `thresholds` sets the purity columns — any number of them, not just two.
    `windows` pins the extent shown per site, as {site: (x0, y0, x1, y1)} in
    LV95; without it the densest patch of parcels is chosen automatically, which
    is convenient but not reproducible between runs if the mask changes.
    `labels` overrides the row captions, e.g. {"rheintal": "St. Galler Rheintal"}.
    """
    import matplotlib.pyplot as plt
    import rasterio
    from matplotlib.patches import Rectangle
    from rasterio.windows import from_bounds

    from .config import CFG, figures_dir, gee_dir

    sites = list(sites or vars(CFG.sites))
    band = list(GROUPS).index(group) + 1
    npx = max(int(zoom_m / RES), 20)

    if master is None:
        print("  reading the parcel set for the polygon column...")
        master = master_parcels(verbose=False)
    codes = GROUPS[group]
    parcels = master[master["_code"].isin(codes)]

    # Locate the imagery before choosing windows. If a scene covers only part
    # of a site, the densest-parcel window can land outside it and the panel
    # comes back black — so the search is restricted to where imagery exists.
    from .config import ROOT
    scenes = {}
    for s in sites:
        p = (rgb or {}).get(s)
        if p is None:
            # Only files that announce themselves as imagery. An earlier
            # version matched any TIFF containing the site name and quietly
            # picked up a rendered z-score map instead.
            cand = (sorted(gee_dir(s).glob("rgb*.tif"))
                    + sorted((ROOT / CFG.paths.gee).glob(f"rgb*{s}*.tif"))
                    + sorted((ROOT / CFG.paths.gee).glob("rgb*.tif")))
            cand = [c for c in cand if s in c.stem or len(sites) == 1]
            p = cand[0] if cand else None
        if p is None:
            continue
        try:
            with rasterio.open(p) as d2:
                scenes[s] = (p, d2.bounds, d2.count, d2.dtypes[0])
                if verbose:
                    a = d2.read(1, out_shape=(1, min(256, d2.height),
                                              min(256, d2.width)))
                    a = a[np.isfinite(a)]
                    print(f"  {s}: {Path(p).name}, {d2.count} bands, "
                          f"{d2.dtypes[0]}, band 1 range "
                          f"{np.nanmin(a):.4g}..{np.nanmax(a):.4g}")
        except Exception as exc:
            print(f"  [!] {s}: cannot open {p} ({exc})")

    rows = []
    with rasterio.open(mask_path()) as ds:
        rb = ds.bounds
        for s in sites:
            try:
                sx0, sy0, sx1, sy1 = _site_bounds_lv95(s)
            except Exception:
                print(f"  [!] unknown site {s!r}"); continue
            x0, y0 = max(rb.left, sx0), max(rb.bottom, sy0)
            x1, y1 = min(rb.right, sx1), min(rb.top, sy1)
            if s in scenes and not (windows or {}).get(s):
                sb = scenes[s][1]
                ix0, iy0 = max(x0, sb.left), max(y0, sb.bottom)
                ix1, iy1 = min(x1, sb.right), min(y1, sb.top)
                if ix1 - ix0 > zoom_m and iy1 - iy0 > zoom_m:
                    x0, y0, x1, y1 = ix0, iy0, ix1, iy1
                else:
                    print(f"  [!] {s}: the scene overlaps the site by less "
                          "than one window; using the full site extent")
            if x1 <= x0 or y1 <= y0:
                print(f"  [!] {s} lies outside the raster"); continue
            pin = (windows or {}).get(s)
            if pin:
                # An explicit window: clipped to the raster so a typo cannot
                # silently read outside it.
                px0, py0 = max(rb.left, pin[0]), max(rb.bottom, pin[1])
                px1, py1 = min(rb.right, pin[2]), min(rb.top, pin[3])
                if px1 <= px0 or py1 <= py0:
                    print(f"  [!] {s}: the given window misses the raster")
                    continue
                zoom = ds.read(band, window=from_bounds(px0, py0, px1, py1,
                                                        ds.transform))
                wx0, wy1 = px0, py1
            else:
                big = ds.read(band,
                              window=from_bounds(x0, y0, x1, y1, ds.transform))
                # densest patch, so the panel shows fields rather than empty
                # ground. Convenient, but it moves if the mask changes — pin it
                # with `windows` for a figure you need to reproduce.
                occ = (big > 0).astype(float)
                best, br, bc = -1, 0, 0
                step = max(npx // 3, 1)
                for r_ in range(0, max(big.shape[0] - npx, 1), step):
                    for c_ in range(0, max(big.shape[1] - npx, 1), step):
                        v = occ[r_:r_ + npx, c_:c_ + npx].sum()
                        if v > best:
                            best, br, bc = v, r_, c_
                zoom = big[br:br + npx, bc:bc + npx]
                wx0 = x0 + bc * RES
                wy1 = y1 - br * RES
            ext = (wx0, wx0 + zoom.shape[1] * RES,
                   wy1 - zoom.shape[0] * RES, wy1)
            rows.append({"site": s, "purity": zoom, "ext": ext})

    if not rows:
        print("nothing to draw")
        return None

    # Crisp polygon rasterisation, no PSF — the cadastre as drawn. The same
    # parcels also give the reference area for the "kept" figures: clipped to
    # the window, so the denominator is the mapped hectares actually on screen.
    # This is the definition area_check() uses, so one number means one thing
    # across the whole analysis.
    from shapely.geometry import box as _box
    cell_ha = (RES * RES) / 1e4
    for r in rows:
        e = r["ext"]
        shape = r["purity"].shape
        sub = parcels.cx[e[0]:e[1], e[2]:e[3]]
        r["poly"] = (_burn(sub.geometry, (e[0], e[3], RES), shape)
                     if len(sub) else np.zeros(shape, np.uint8))
        if len(sub):
            clip = _box(e[0], e[2], e[1], e[3])
            keep_ = sub[sub.geometry.intersects(clip)].copy()
            keep_["geometry"] = keep_.geometry.intersection(clip)
            keep_ = keep_[keep_.geometry.area > 0]
            r["parcels"] = keep_
            r["poly_ha"] = float(keep_.geometry.area.sum() / 1e4)
        else:
            r["parcels"] = sub
            r["poly_ha"] = 0.0

    # imagery, if it was exported
    for r in rows:
        r["rgb"] = None
        if r["site"] not in scenes:
            continue
        p, _b, ncount, _dt = scenes[r["site"]]
        idx = list(rgb_bands)
        if ncount < max(idx):
            print(f"  [!] {r['site']}: the file has {ncount} band(s) but "
                  f"bands {idx} were requested — export B8/B4/B3 as three "
                  "bands, or set rgb_bands")
            continue
        try:
            with rasterio.open(p) as ds:
                e = r["ext"]
                w = from_bounds(e[0], e[2], e[1], e[3], ds.transform)
                a = ds.read(indexes=idx, window=w,
                            out_shape=(3, r["purity"].shape[0],
                                       r["purity"].shape[1]),
                            boundless=True, fill_value=0).astype(np.float32)
            img = np.moveaxis(a, 0, -1)
            # Percentile stretch PER BAND on its own valid pixels. Sentinel-2
            # L2A comes as integer reflectance x10000, so a shared stretch
            # across bands would leave the image looking like a grey ramp —
            # NIR is far brighter than red or green over vegetation, and it is
            # exactly that difference the false colour is meant to show.
            for b in range(img.shape[-1]):
                ch = img[..., b]
                finite = ch[np.isfinite(ch) & (ch > 0)]
                if finite.size < 16:
                    img[..., b] = 0.0
                    continue
                lo, hi = np.quantile(finite, stretch)
                if hi - lo < 1e-9:
                    lo, hi = float(finite.min()), float(finite.max())
                img[..., b] = np.clip((ch - lo) / max(hi - lo, 1e-9), 0, 1)
            if verbose:
                mn = [float(np.nanmean(img[..., b])) for b in range(img.shape[-1])]
                print(f"  {r['site']}: stretched band means "
                      + ", ".join(f"{v:.2f}" for v in mn))
            nonzero = float((img.sum(axis=-1) > 0.01).mean())
            if nonzero < 0.5:
                print(f"  [!] {r['site']}: the exported scene covers only "
                      f"{100 * nonzero:.0f}% of the window being drawn — "
                      "re-export over the whole site, or pass a smaller "
                      "zoom_m")
            r["rgb"] = img
        except Exception as exc:
            print(f"  [!] {r['site']}: could not read {p} ({exc})")

    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Polygon as MplPoly

    plt.style.use("dark_background")
    nr, nc = len(rows), 2 + len(thresholds)
    fig, axes = plt.subplots(nr, nc, figsize=(2.9 * nc, 3.0 * nr + 0.9),
                             dpi=180, squeeze=False)
    heads = ["Sentinel-2, false colour", "Parcel polygons"] + \
            [f"Purity \u2265 {t:g}" for t in thresholds]

    for i, r in enumerate(rows):
        panels = [r["rgb"], r["poly"]] + \
                 [np.where(r["purity"] >= max(t, 1e-9), r["purity"], np.nan)
                  for t in thresholds]
        for j, dat in enumerate(panels):
            ax = axes[i][j]
            e = r["ext"]
            if j == 0:
                if dat is None:
                    ax.set_facecolor("#0d1117")
                    ax.text(0.5, 0.5, "no imagery exported\nfor this site",
                            ha="center", va="center", fontsize=8.5,
                            color="#6e7681", transform=ax.transAxes)
                else:
                    ax.imshow(dat, interpolation="nearest", extent=e)
            elif j == 1:
                # VECTOR polygons, not a raster. The point of this panel is the
                # sharpness of the boundary — rasterising it here would show the
                # same stair-stepping the purity panels have, and lose the
                # contrast the column exists to make.
                ax.set_facecolor("#000000")
                pats = []
                for geom in r["parcels"].geometry:
                    if geom is None or geom.is_empty:
                        continue
                    gs = geom.geoms if geom.geom_type == "MultiPolygon" else [geom]
                    for gg in gs:
                        pats.append(MplPoly(np.asarray(gg.exterior.coords),
                                            closed=True))
                if pats:
                    ax.add_collection(PatchCollection(
                        pats, facecolor="white", edgecolor="#9aa0a6",
                        linewidths=0.35))
                ax.set_xlim(e[0], e[1]); ax.set_ylim(e[2], e[3])
                ax.set_aspect("equal")
            else:
                ax.imshow(np.ma.masked_invalid(dat), cmap="gray", vmin=0,
                          vmax=1, interpolation="nearest", extent=e)
                keep = r["purity"] >= max(thresholds[j - 2], 1e-9)
                pct = (100 * float(keep.sum()) * cell_ha / r["poly_ha"]
                       if r["poly_ha"] > 0 else float("nan"))
                red = abs(thresholds[j - 2] - PURITY_MIN) < 1e-9
                ax.text(0.97, 0.965, f"{pct:.0f}% of mapped area",
                        transform=ax.transAxes,
                        ha="right", va="top", fontsize=9, fontweight="bold",
                        color="#ff3b30" if red else "#c9d1d9",
                        bbox=dict(facecolor="black", alpha=0.72,
                                  edgecolor="none", boxstyle="round,pad=0.25"))

            ax.set_xticks([]); ax.set_yticks([])
            if i == 0:
                ax.set_title(heads[j], fontsize=10, fontweight="bold",
                             color="white", loc="left", pad=6)
            if j == 0:
                nm = (labels or {}).get(r["site"])
                if nm is None:
                    sc = getattr(CFG.sites, r["site"], None)
                    nm = getattr(sc, "label", None) or r["site"].capitalize()
                ax.set_ylabel(nm, fontsize=11, fontweight="bold",
                              color="#ff3b30", labelpad=8)
            for sp in ax.spines.values():
                sp.set_color("#30363d")
        # scale bar on the last panel of each row
        e = r["ext"]
        w_ = e[1] - e[0]
        bx, by = e[0] + 0.05 * w_, e[2] + 0.05 * (e[3] - e[2])
        axL = axes[i][-1]
        axL.add_patch(Rectangle((bx - 0.012 * w_, by - 0.012 * w_),
                                100 + 0.024 * w_, 0.062 * w_,
                                color="black", alpha=0.72, zorder=5))
        axL.add_patch(Rectangle((bx, by), 100, 0.009 * w_, color="white",
                                zorder=6))
        axL.text(bx + 50, by + 0.016 * w_, "100 m", ha="center", va="bottom",
                 fontsize=7, color="white", zorder=6)

    fig.suptitle("From imagery to mask: what each step keeps",
                 fontsize=12.5, fontweight="bold", color="#ff3b30",
                 x=0.012, ha="left", y=0.985)
    fig.supxlabel(f"{group} grassland at {RES} m, {zoom_m} m across. Purity is "
                  "the share of a cell's sensor response coming from inside a "
                  f"parcel ({psf_caption()}). White "
                  "is pure, grey is mixed, black is discarded.\n"
                  "Percentages are retained area against the polygon area in "
                  "the same window \u2014 above 100% at a low threshold means "
                  "the response is reaching past the parcels.",
                  fontsize=7.5, color="#6e7681")
    fig.tight_layout(rect=[0.004, 0.055, 0.998, 0.958])
    fig.subplots_adjust(hspace=0.06, wspace=0.035)
    out = out or (figures_dir("mask") / f"mask_overview_{group}.png")
    fig.savefig(out, dpi=180, facecolor=fig.get_facecolor())
    print(f"  -> {out}")
    plt.show()
    return rows


def figure_parcel_size(master=None, group=None, out=None, year=None,
                       bins=(0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0), shape_clip=16):
    """How big and how convoluted the parcels are — the case for slide A.

    Two panels, because size alone does not explain what the mask costs.

    LEFT   parcel area, log scale, one violin per land-use group. The number a
           lay audience needs: a 10 m pixel is 0.01 ha, so a median field holds
           about a hundred of them. That sounds comfortable.
    RIGHT  the shape index, perimeter over the square root of area. 4.0 is a
           perfect square; Swiss grassland runs 5 to 6.5. That is where the
           comfort goes: a long narrow field is mostly edge, however many
           pixels fit inside it.

    Printed alongside: the share of parcels under half a hectare, for which a
    hundred pixels becomes a dozen.
    """
    import matplotlib.pyplot as plt

    from .config import figures_dir

    if master is None:
        print("  building the parcel set — pass master= to skip this")
        master = master_parcels(verbose=False)
    groups = [group] if group else ["intensive", "extensive", "mown", "grazed"]

    data = {}
    for n in groups:
        sub = master[master["_code"].isin(GROUPS[n])]
        if len(sub) < 10:
            continue
        a = (sub.geometry.area / 1e4).to_numpy()
        shp = (sub.geometry.length / np.sqrt(sub.geometry.area)).to_numpy()
        ok = np.isfinite(a) & np.isfinite(shp) & (a > 0)
        data[n] = (a[ok], shp[ok])
    if not data:
        print("no group has enough parcels")
        return None
    names = list(data)
    col = {"intensive": "#3fb950", "extensive": "#58a6ff",
           "mown": "#d29922", "grazed": "#8250df"}

    plt.style.use("dark_background")
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.2), dpi=180,
                             constrained_layout=True)

    ax = axes[0]
    vals = [np.log10(data[n][0]) for n in names]
    parts = ax.violinplot(vals, positions=range(len(names)), widths=0.85,
                          showextrema=False, showmedians=True)
    for b, n, v in zip(parts["bodies"], names, vals):
        b.set_facecolor(col[n]); b.set_alpha(0.78); b.set_edgecolor("none")
        vt = b.get_paths()[0].vertices
        vt[:, 1] = np.clip(vt[:, 1], v.min(), v.max())
    parts["cmedians"].set_color("white"); parts["cmedians"].set_linewidth(1.9)
    ax.set_yticks([np.log10(b) for b in bins])
    ax.set_yticklabels([f"{b:g}" for b in bins])
    ax.set_ylabel("parcel area [ha], log scale", fontsize=10, color="white")
    ax.set_title("How big is a field?", fontsize=11.5, fontweight="bold",
                 color="white", loc="left")
    # what a 10 m pixel is, on the same axis
    ax.axhline(np.log10(0.01), color="#eda100", ls="--", lw=1.5)
    ax.annotate("one 10 m pixel = 0.01 ha", xy=(0.985, np.log10(0.01)),
                xycoords=("axes fraction", "data"), xytext=(0, 6),
                textcoords="offset points", ha="right", fontsize=9,
                color="#eda100")

    ax = axes[1]
    # The shape index has a long upper tail — a minority of parcels run well
    # past 10 — so the violins are clipped for display. Clipping too low makes
    # every violin end in a flat top, which reads as data at the ceiling rather
    # than as a cut-off; the axis says so when the clip actually bites.
    n_clipped = int(sum(int((data[n][1] > shape_clip).sum()) for n in names))
    vals2 = [np.clip(data[n][1], 3.6, shape_clip) for n in names]
    parts = ax.violinplot(vals2, positions=range(len(names)), widths=0.85,
                          showextrema=False, showmedians=True)
    for b, n, v in zip(parts["bodies"], names, vals2):
        b.set_facecolor(col[n]); b.set_alpha(0.78); b.set_edgecolor("none")
        vt = b.get_paths()[0].vertices
        vt[:, 1] = np.clip(vt[:, 1], v.min(), v.max())
    parts["cmedians"].set_color("white"); parts["cmedians"].set_linewidth(1.9)
    ax.axhline(4.0, color="#eda100", ls="--", lw=1.5)
    ax.annotate("4.0 = a perfect square", xy=(0.985, 4.0),
                xycoords=("axes fraction", "data"), xytext=(0, 6),
                textcoords="offset points", ha="right", fontsize=9,
                color="#eda100")
    ax.set_ylabel("perimeter / \u221aarea" +
                  (f"  (display clipped at {shape_clip:g})" if n_clipped else ""),
                  fontsize=10, color="white")
    ax.set_title("How much of it is edge?", fontsize=11.5, fontweight="bold",
                 color="white", loc="left")

    for ax_ in axes:
        ax_.set_xticks(range(len(names)))
        ax_.set_xticklabels(
            [f"{n}\nn={len(data[n][0]):,}" for n in names], fontsize=9)
        ax_.grid(alpha=0.22, ls=":", color="gray")
        ax_.tick_params(colors="white", labelsize=9)
        for sp in ("top", "right"):
            ax_.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax_.spines[sp].set_color("#484f58")

    print(f"  {'group':<12}{'parcels':>9}{'median ha':>11}{'<0.5 ha':>9}"
          f"{'median shape':>14}{'pixels in a median field':>26}")
    for n in names:
        a, shp = data[n]
        print(f"  {n:<12}{len(a):>9,}{np.median(a):>11.2f}"
              f"{100 * (a < 0.5).mean():>8.0f}%{np.median(shp):>14.2f}"
              f"{np.median(a) / 0.01:>26.0f}")
    if n_clipped:
        tot = sum(len(data[n][0]) for n in names)
        print(f"\n  {n_clipped:,} of {tot:,} parcels ({100 * n_clipped / tot:.1f}%) "
              f"have a shape index above {shape_clip:g} and are clipped in the "
              "right panel only.")

    fig.suptitle("Swiss grassland parcels: small, and mostly edge",
                 fontsize=12.5, fontweight="bold", color="#ff3b30",
                 x=0.012, ha="left")
    fig.supxlabel("A median field holds about a hundred 10 m pixels — but a "
                  "shape index near 6 means much of it lies within a pixel or "
                  "two of a boundary.", fontsize=7.5, color="#6e7681")
    out = out or (figures_dir("mask") / "parcel_size.png")
    fig.savefig(out, dpi=180, facecolor=fig.get_facecolor())
    print(f"\n  -> {out}")
    plt.show()
    return data


def figure_psf(out=None, mtf=None, span=40, ground=True):
    """What one pixel actually sees — cross-section, surface, and the ground.

    The figure the whole purity argument rests on, for an audience that has
    never met a point-spread function.

    LEFT   the response in cross-section, with the nominal pixel drawn to the
           same scale beneath it. The single fact to take away: the response is
           wider than the grid it is sampled on.
    RIGHT  the same response as a surface over `span` metres of ground, with
           the pixel marked on the floor and, below it, a patch of farmland
           crossed by a brook. The pixel covers parts of two fields and clips
           the brook — while the whole area still contributes.

    The response is NOT a bare Gaussian. It is the detector footprint, a box of
    one pixel width, convolved with a Gaussian optical response — which is why
    it has a flatter top than a Gaussian would. Sigma is derived, not assumed:
    the detector's own MTF at Nyquist is sinc(1/2) = 0.637, so to reach a
    specified system MTF of 0.25 the optics must supply 0.393, and that fixes
    sigma at about 4.35 m.

    Approximations worth stating in the caption: the optics are Airy-like
    rather than Gaussian; the real response is slightly wider along-track
    because of motion during integration; and ESA publishes MTF per band, while
    one representative value is used here.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as MplPoly
    from scipy.special import erf

    from .config import figures_dir

    pp = psf_params(mtf=mtf, model="box_gauss")
    m_, mtf_det, mtf_opt = pp["mtf_system"], pp["mtf_detector"], pp["mtf_optics"]
    sig, fwhm = pp["sigma_m"], pp["fwhm_m"]

    def prof(x):
        return 0.5 * (erf((x + RES / 2) / (sig * np.sqrt(2)))
                      - erf((x - RES / 2) / (sig * np.sqrt(2))))

    pk = prof(0.0)
    xs = np.linspace(-span / 2, span / 2, 1200)
    ys = prof(xs) / pk

    plt.style.use("dark_background")
    fig = plt.figure(figsize=(13.4, 6.4), dpi=180)
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.18], wspace=0.10,
                          left=0.055, right=0.985, top=0.87, bottom=0.20)

    # ---- cross-section ---------------------------------------------------
    ax = fig.add_subplot(gs[0])
    ax.axvspan(-RES / 2, RES / 2, color="#eda100", alpha=0.13, lw=0)
    for e in (-RES / 2, RES / 2):
        ax.axvline(e, color="#eda100", ls=":", lw=1.2)
    ax.fill_between(xs, ys, color="white", alpha=0.10, lw=0)
    ax.plot(xs, ys, color="white", lw=2.6)
    ax.plot([-fwhm / 2, fwhm / 2], [0.5, 0.5], color="#eda100", lw=1.8)
    for e in (-fwhm / 2, fwhm / 2):
        ax.plot([e, e], [0.46, 0.54], color="#eda100", lw=1.8)
    ax.annotate(f"FWHM {fwhm:.1f} m", xy=(fwhm / 2, 0.5),
                xytext=(fwhm / 2 + 2.5, 0.60), fontsize=10, color="#eda100",
                arrowprops=dict(arrowstyle="-", color="#eda100", lw=0.9))
    ax.annotate(f"the nominal {RES:g} m pixel", xy=(RES / 2, 0.93),
                xytext=(RES / 2 + 3, 1.02), fontsize=10, color="#eda100",
                arrowprops=dict(arrowstyle="-", color="#eda100", lw=0.9))
    ax.set_xlim(-span / 2, span / 2); ax.set_ylim(0, 1.16)
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_xlabel("metres from the pixel centre", fontsize=10, color="white")
    ax.set_ylabel("response, normalised to peak", fontsize=10, color="white")
    ax.set_title("In cross-section", fontsize=11.5, fontweight="bold",
                 color="white", loc="left")
    ax.grid(alpha=0.18, ls=":", color="gray")

    # ---- surface over ground --------------------------------------------
    ax2 = fig.add_subplot(gs[1])
    ax2.set_axis_off()
    yaw, pit, S = -0.70, 0.56, 1.0
    cy, sy = np.cos(yaw), np.sin(yaw)
    cp, sp = np.cos(pit), np.sin(pit)

    def P(u, v, w):
        # NOTE the sign on the vertical: matplotlib's y axis points up while a
        # canvas points down. Porting the projection without flipping it put
        # the surface in a pit and the ground plane in the sky.
        x, y, z = u - 0.5, v - 0.5, w - 0.5
        return ((x * cy - y * sy) * 0.82 * S,
                -((x * sy + y * cy) * cp - z * sp) * 0.42 * S)

    def draw(pts, fc=None, ec=None, lw=1.0, z=1):
        ax2.add_patch(MplPoly(np.array(pts), closed=True,
                              facecolor=fc or "none", edgecolor=ec or "none",
                              lw=lw, zorder=z))

    GZ = -1.95 if ground else 0.0
    if ground:
        gq = [P(0, 0, GZ), P(1, 0, GZ), P(1, 1, GZ), P(0, 1, GZ)]
        draw(gq, fc="#3d4a32", z=1)
        bnd = [0.0, 0.26, 0.50, 0.74, 1.0]
        tint = ["#6f8348", "#93a05c", "#b0a068", "#5c7340"]
        for k in range(4):
            draw([P(0, bnd[k], GZ), P(1, bnd[k] - 0.08, GZ),
                  P(1, bnd[k + 1] - 0.08, GZ), P(0, bnd[k + 1], GZ)],
                 fc=tint[k], z=2)
        uu = np.linspace(0, 1, 12)
        rv = 0.66 - 0.26 * uu
        riv = [P(u, r + 0.028, GZ) for u, r in zip(uu, rv)] + \
              [P(u, r - 0.028, GZ) for u, r in zip(uu[::-1], rv[::-1])]
        draw(riv, fc="#46687c", z=3)
        draw(gq, ec="#ffffff55", lw=1.2, z=4)

    hp = (RES / span) / 2
    sq = [(0.5 - hp, 0.5 - hp), (0.5 + hp, 0.5 - hp),
          (0.5 + hp, 0.5 + hp), (0.5 - hp, 0.5 + hp)]
    if ground:
        draw([P(u, v, GZ) for u, v in sq], fc="#eda10033", ec="#eda100",
             lw=1.8, z=5)
        for u, v in sq:
            a, b = P(u, v, 0), P(u, v, GZ)
            ax2.plot([a[0], b[0]], [a[1], b[1]], color="#eda100", lw=0.9,
                     ls=(0, (3, 3)), zorder=5)

    draw([P(0, 0, 0), P(1, 0, 0), P(1, 1, 0), P(0, 1, 0)],
         fc="#ffffff08", ec="#ffffff6b", lw=1.3, z=6)
    for t in (0.25, 0.5, 0.75):
        for a, b in ((P(t, 0, 0), P(t, 1, 0)), (P(0, t, 0), P(1, t, 0))):
            ax2.plot([a[0], b[0]], [a[1], b[1]], color="#ffffff33", lw=0.9,
                     zorder=6)

    N = 30
    g1 = np.linspace(-span / 2, span / 2, N + 1)
    f1 = prof(g1) / pk
    quads = []
    for i in range(N):
        for j in range(N):
            pts = [P(i / N, j / N, f1[i] * f1[j]),
                   P((i + 1) / N, j / N, f1[i + 1] * f1[j]),
                   P((i + 1) / N, (j + 1) / N, f1[i + 1] * f1[j + 1]),
                   P(i / N, (j + 1) / N, f1[i] * f1[j + 1])]
            d = ((i / N - 0.5) * sy + (j / N - 0.5) * cy) * sp \
                + (f1[i] * f1[j] - 0.5) * cp
            quads.append((d, pts))
    for _d, pts in sorted(quads, key=lambda q: q[0]):
        draw(pts, fc="#000000", ec="#ffffffd0", lw=0.45, z=7)

    draw([P(u, v, 0) for u, v in sq], fc="#eda10055", ec="#eda100", lw=2,
         z=8)
    c = P(0.5 + hp, 0.5 - hp, 0)
    ax2.annotate(f"{RES:g} \u00d7 {RES:g} m pixel", xy=c,
                 xytext=(c[0] + 0.12, c[1] - 0.10), fontsize=9.5,
                 color="#eda100",
                 arrowprops=dict(arrowstyle="-", color="#eda100", lw=0.9))
    # Edge labels offset outward from the floor centre, so they never land on
    # the mesh.
    fc = P(0.5, 0.5, 0)
    for (a_, b_) in ((P(0, 0, 0), P(1, 0, 0)), (P(1, 0, 0), P(1, 1, 0))):
        mx, my = (a_[0] + b_[0]) / 2, (a_[1] + b_[1]) / 2
        dx, dy = mx - fc[0], my - fc[1]
        n = max((dx * dx + dy * dy) ** 0.5, 1e-6)
        ax2.plot([a_[0] + dx / n * 0.05, b_[0] + dx / n * 0.05],
                 [a_[1] + dy / n * 0.05, b_[1] + dy / n * 0.05],
                 color="#ffffff70", lw=1, zorder=9)
        ax2.text(mx + dx / n * 0.12, my + dy / n * 0.12, f"{span:g} m",
                 fontsize=9, color="#ffffffb0", ha="center", va="center",
                 zorder=9)
    if ground:
        gl = P(0, 0.5, GZ)
        ax2.text(gl[0] - 0.03, gl[1], "the ground\nbelow", fontsize=9.5,
                 color="#ffffffcc", ha="right", va="center")
    ax2.set_title(f"Over {span:g} \u00d7 {span:g} m of ground",
                  fontsize=11.5, fontweight="bold", color="white", loc="left",
                  pad=12)
    # Limits fitted to the content, so the panel is not shrunk by empty space.
    ax2.set_xlim(-0.66, 0.78)
    ax2.set_ylim(-0.86 if ground else -0.30, 0.30)
    ax2.set_aspect("equal")

    for a_ in (ax,):
        a_.tick_params(colors="white", labelsize=9)
        for sp_ in ("top", "right"):
            a_.spines[sp_].set_visible(False)
        for sp_ in ("left", "bottom"):
            a_.spines[sp_].set_color("#484f58")

    vals = [("\u03c3 optics", f"{sig:.2f} m"), ("FWHM", f"{fwhm:.1f} m"),
            ("detector MTF", f"{mtf_det:.3f}"),
            ("optics MTF", f"{mtf_opt:.3f}"), ("system MTF", f"{m_:.2f}")]
    for i, (k, v) in enumerate(vals):
        x = 0.055 + i * 0.185
        fig.text(x, 0.115, k, fontsize=8.5, color="#8b949e")
        fig.text(x, 0.072, v, fontsize=15, fontweight="bold", color="white")

    fig.suptitle("A pixel is blurrier than its own grid", fontsize=13,
                 fontweight="bold", color="#ff3b30", x=0.012, ha="left",
                 y=0.965)
    fig.text(0.055, 0.018,
             "MTF at Nyquist. The detector footprint, a box one pixel wide, "
             "convolved with a Gaussian optical response \u2014 which is why "
             "the top is flatter than a Gaussian.\nApproximations: the optics "
             "are Airy-like rather than Gaussian, the real response is "
             "slightly wider along-track, and one representative MTF stands "
             "for all bands.", fontsize=7, color="#6e7681")
    out = out or (figures_dir("mask") / "psf.png")
    fig.savefig(out, dpi=180, facecolor=fig.get_facecolor())
    print(f"  sigma {sig:.2f} m, FWHM {fwhm:.1f} m  ->  {out}")
    plt.show()
    return {"sigma": float(sig), "fwhm": float(fwhm),
            "mtf_detector": float(mtf_det), "mtf_optics": float(mtf_opt)}
