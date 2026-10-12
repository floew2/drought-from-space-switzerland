"""Which fields: reading the annual deliveries and building the master set.

Parcels are matched between years by AREA OVERLAP, never by identity: outlines
are redrawn, split and merged between deliveries, and IDs are not stable.
Persistence selects parcels and is never applied to a raster. METHOD §2.1.
"""

from __future__ import annotations

import re
from pathlib import Path

import geopandas as gpd
import pandas as pd

from ..config import parcels_dir
from .settings import (
    CODE_FIELD_CANDIDATES,
    CRS,
    GRANDFATHER_CODES,
    GROUPS,
    MIN_OVERLAP,
    MIN_YEARS_SAME,
    PATTERN,
    SIMPLIFY_M,
    SLIVER_FLOOR_M2,
)

__all__ = ["category_history", "category_of", "load_year", "master_parcels",
           "matching_check"]


def _year_files(shp_dir=None) -> dict:
    d = Path(shp_dir) if shp_dir else parcels_dir()
    out = {}
    for f in sorted(d.rglob(PATTERN)):
        m = re.search(r"(20\d\d)", f.stem)
        if m:
            out.setdefault(int(m.group(1)), f)
    if not out:
        raise SystemExit(
            f"No files matching {PATTERN} under {d.resolve()} (searched "
            "recursively).\nPut the yearly Nutzungsflaechen shapefiles there — "
            "one folder per year is fine — or pass shp_dir.")
    return out


def load_year(path: Path) -> gpd.GeoDataFrame:
    """Read one annual delivery: valid polygons with an integer `_code`.

    The CRS is asserted, not assumed: a file without one is refused rather than
    labelled LV95 on faith. A file in another CRS is reprojected.

    Args:
        path: Shapefile of one delivery.

    Returns:
        GeoDataFrame in EPSG:2056 with the land-use code as `_code`.

    Raises:
        ValueError: If the file carries no CRS.
        SystemExit: If no land-use code field is found.
    """
    g = gpd.read_file(path)
    if g.crs is None:
        raise ValueError(f"{path.name} has no CRS; expected {CRS}. Add the .prj "
                         "rather than assuming one.")
    if g.crs.to_epsg() != 2056:
        print(f"    reprojecting {path.name} from {g.crs.to_string()} to {CRS}")
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
    """Make every geometry valid without losing area, or drop it with a count.

    make_valid first, buffer(0) only as a fallback. The order matters:
    buffer(0) "repairs" a figure-eight outline into a valid polygon with one
    lobe missing, half the field gone without any error, so trying it first
    meant make_valid never saw that case. make_valid keeps every lobe. A ring
    that touches itself at a point is also handled; left alone, it throws a
    side-location conflict several minutes into the overlay. Anything still
    broken is discarded with a count: one unusable parcel is not worth losing
    the overlay for.
    """
    from shapely.validation import make_valid

    bad = ~g.geometry.is_valid
    n = int(bad.sum())
    if n:
        fixed = g.loc[bad, "geometry"].apply(make_valid)
        still = ~fixed.is_valid
        if still.any():
            fixed.loc[still] = fixed.loc[still].buffer(0)
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


def category_of(code: object) -> str | None:
    """Management category of a land-use code; None if it is not grassland.

    intensive and extensive are checked first: they are the analytical split.
    mown / grazed is the broader "is there a cut signal" split, and the two do
    not nest, so a code can belong to both.

    Args:
        code: LNF land-use code (anything int() accepts).

    Returns:
        "intensive", "extensive", "grazed", "mown_other" or None.
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


def master_parcels(shp_dir: str | Path | None = None, verbose: bool = True,
                   min_overlap: float | None = None, min_years: int | None = None,
                   sliver_floor_m2: float | None = None) -> gpd.GeoDataFrame:
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
        The master parcel set, EPSG:2056, with `_code` and `_cat`. The funnel
        (counts dropped and kept, per category) is in `.attrs["funnel"]`.
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
    lost_grass = int(((n_grass < min_years) & ~grandfathered).sum())
    lost_cat = int(((n_grass >= min_years) & (n_same < min_years) & ~grandfathered).sum())
    held_all = int(((n_same >= len(years)) & ~grandfathered).sum())
    out.attrs["funnel"] = {
        "year": int(cur_year), "grassland": int(n0), "sliver_floor_m2": float(floor),
        "dropped_not_grassland": lost_grass, "dropped_changed_category": lost_cat,
        "exempt": int(grandfathered.sum()), "remain": len(out),
        "remain_share": round(len(out) / max(n0, 1), 4),
        "remain_ha": round(float(out.geometry.area.sum() / 1e4), 1),
        "held_category_all_years_share": round(held_all / max(n0, 1), 4),
        "per_category": {n: {"parcels": int(out["_code"].isin(c).sum()),
                             "ha": round(float(out[out["_code"].isin(c)].geometry.area.sum()
                                               / 1e4), 1)}
                         for n, c in GROUPS.items()},
    }
    if verbose:
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


def category_history(shp_dir: str | Path | None = None, min_overlap: float | None = None,
                     verbose: bool = True, sliver_floor_m2: float | None = None) -> pd.DataFrame:
    """How many years each current parcel held its category, by area overlap.

    The companion to `master_parcels`: the same matching, reporting the whole
    distribution instead of applying a cut-off, so the rule can be chosen from
    evidence. It counts separately how often a parcel was ANY grassland and how
    often it was ITS OWN category; the gap between the two is management change.

    Args:
        shp_dir: Folder with the yearly shapefiles (default: config).
        min_overlap: Share of a field that must lie inside a year's category.
        verbose: Print the distribution and what each rule would keep.
        sliver_floor_m2: Override the sliver floor (2500 = old behaviour).

    Returns:
        DataFrame per current parcel: area_ha, cat, years_grassland,
        years_same_cat, grandfathered.
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


def matching_check(shp_dir: str | Path | None = None,
                   overlaps: tuple[float, ...] = (0.5, 0.7, 0.8, 0.9, 0.95)) -> dict:
    """How sensitive the master set is to the overlap threshold.

    A gentle slope means parcels are either clearly the same ground or clearly
    not, and the matching is robust; a steep one means the threshold is doing
    more work than it should.

    Args:
        shp_dir: Folder with the yearly shapefiles (default: config).
        overlaps: Overlap thresholds to test.

    Returns:
        dict threshold -> % of parcels that held their category in all years.
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
