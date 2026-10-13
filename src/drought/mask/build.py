"""Building the mask, and the numbers that check it. METHOD §2.4.

`build_mask` writes the purity raster; `area_table` measures what each purity
threshold costs; `psf_model_check` measures whether the response model changes
that. The ratio at purity > 0 must exceed 1 for every category.
"""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

from ..config import mask_path
from .io import site_bounds_lv95, write_raster
from .parcels import master_parcels
from .purity import psf_caption, psf_params, purity_psf
from .settings import CLIP_BBOX_LV95, CRS, GROUPS, PURITY_MIN, RES, SUBPIXEL

__all__ = ["area_table", "build_mask", "main", "output_grid", "psf_model_check",
           "window_polygon_ha",
           "verify_mask"]


def output_grid(master: gpd.GeoDataFrame, verbose: bool = True) -> tuple[tuple, tuple]:
    """The 10 m output grid for a parcel set, padded so the halo is never cut.

    Args:
        master: Parcel set (GeoDataFrame in EPSG:2056).
        verbose: Say when the grid is clipped to `mask.clip_bbox_lv95`.

    Returns:
        ((minx, maxy, res), (height, width)).
    """
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


def window_polygon_ha(geoms, grid: tuple, shape: tuple) -> float:
    """Polygon area [ha] inside the raster window.

    The reference for an area ratio must cover the same ground as the raster.
    With `mask.clip_bbox_lv95` set, parcels reach beyond the grid; their area
    outside it exists as polygon but can never appear as pixels (bug 13).

    Args:
        geoms: Parcel geometries (EPSG:2056).
        grid: (minx, maxy, res) as returned by `output_grid`.
        shape: (height, width) of the raster.

    Returns:
        Area in hectares.
    """
    from shapely.geometry import box as _box

    minx, maxy, res = grid
    win = _box(minx, maxy - shape[0] * res, minx + shape[1] * res, maxy)
    g = gpd.GeoSeries(geoms).intersection(win)
    return float(g.area.sum() / 1e4)


def psf_model_check(master: gpd.GeoDataFrame, thresholds: tuple[float, ...] = (0.0, 0.5, 0.9),
                    verbose: bool = True) -> pd.DataFrame:
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


def build_mask(shp_dir: str | Path | None = None,
               out_tif: str | Path | None = None) -> gpd.GeoDataFrame:
    """Build the four-band purity mask from the yearly deliveries.

    Master parcel set, then PSF purity per category on the 10 m grid, written
    as float32 purity 0-1. No persistence band: persistence is already in the
    parcel selection. The threshold is applied downstream (in Earth Engine).

    Args:
        shp_dir: Folder with the yearly shapefiles (default: config).
        out_tif: Output path (default: config `paths.masks`).

    Returns:
        The master parcel set the mask was built from.
    """
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
        poly_ha = window_polygon_ha(sub.geometry, grid, (height, width))
        h0 = float((p > 0).sum()) * cell_ha
        h9 = float((p >= PURITY_MIN).sum()) * cell_ha
        print(f"  {name:<12}{len(sub):>9,}{poly_ha:>12,.0f}{h0:>11,.0f}"
              f"{h0 / poly_ha:>7.2f}{h9:>13,.0f}{h9 / poly_ha:>7.2f}")
        bands.append((p, name))

    print("\n  polygon ha: inside the raster window. The ratio at >0 MUST exceed 1:\n"
          "  the PSF response spreads beyond every boundary. If it does not,\n"
          "  something is wrong.")

    out_tif = Path(out_tif) if out_tif else mask_path()
    write_raster(bands, grid, (height, width), out_tif)
    from ..provenance import record, write_sidecar
    from .parcels import _year_files
    deliveries = {str(y): _sha_short(f) for y, f in _year_files(shp_dir).items()}
    write_sidecar(out_tif, record("mask", out_tif, {
        "deliveries_sha256": deliveries, "bands": [n for _, n in bands],
        "parcels": {n: int(master["_code"].isin(c).sum()) for n, c in GROUPS.items()},
        "master_parcels": len(master)}))
    print(f"\n  wrote {out_tif} ({out_tif.stat().st_size / 1e6:.1f} MB)")
    print("  band order: " + ", ".join(f"b{i} {n}"
                                       for i, (_, n) in enumerate(bands, 1)))
    print("\n  FLOAT purity 0-1. No persistence multiplication — persistence is "
          "already in\n  the parcel selection. Threshold in GEE with .gte(0.9), "
          "and upload with\n  pyramidingPolicy {'.default': 'sample'} or reduce "
          "at the native 10 m.")
    return master


main = build_mask        # the name the notebook and older docs use


def _sha_short(path) -> str:
    import hashlib
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def verify_mask() -> dict:
    """Refuse a mask built under other settings than the current ones.

    Returns:
        The mask's provenance record.

    Raises:
        drought.provenance.ProvenanceError: Naming the setting that differs.
    """
    from ..provenance import verify
    return verify(mask_path(), "mask")


def area_table(master: gpd.GeoDataFrame | None = None, site: str | None = None,
               thresholds: tuple[float, ...] = (0.0, 0.25, 0.5, 0.75, 0.9),
               verbose: bool = True) -> dict:
    """Observed area against polygon area, per category and purity threshold.

    Both sides describe the same fixed parcel set, clipped to the ground the
    raster covers (parcels outside the grid exist as polygons but cannot exist
    in the raster). The ratio at > 0 must exceed 1, because the response
    reaches past every boundary; crossing 1 is not a target.

    Args:
        master: Master parcel set; built if None (slow).
        site: Restrict to a site window, or None for the full raster.
        thresholds: Purity thresholds; 0 means "> 0".
        verbose: Print the tables.

    Returns:
        dict category -> {parcels, polygon_ha, median_ha, median_shape,
        observed_ha {t: ha}, ratio {t: ratio}}, plus "_window" (the raster
        window, for plotting).
    """
    import rasterio
    from shapely.geometry import box as _box

    if master is None:
        master = master_parcels(verbose=False)
    verify_mask()
    with rasterio.open(mask_path()) as ds:
        rb = ds.bounds
        if site:
            sx0, sy0, sx1, sy1 = site_bounds_lv95(site)
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

    out = {}
    for n in names:
        sub = master[master["_code"].isin(GROUPS[n])]
        a_ha = sub.geometry.area / 1e4
        shp = (sub.geometry.length / np.sqrt(sub.geometry.area)).median() \
            if len(sub) else float("nan")
        obs = {t: float((arrs[n] >= max(t, 1e-9)).sum()) * cell_ha for t in thresholds}
        poly = float(a_ha.sum())
        out[n] = {"parcels": len(sub), "polygon_ha": poly,
                  "median_ha": float(a_ha.median()) if len(sub) else np.nan,
                  "median_shape": float(shp), "observed_ha": obs,
                  "ratio": {t: (v / poly if poly else np.nan) for t, v in obs.items()}}
    if verbose:
        print(f"{site or 'full extent'}  ·  master parcel set  ·  reference "
              "clipped to the raster extent\n")
        print(f"  {'category':<12}{'parcels':>9}{'polygon ha':>12}{'median ha':>11}"
              f"{'shape':>8}")
        for n in names:
            d = out[n]
            print(f"  {n:<12}{d['parcels']:>9,}{d['polygon_ha']:>12,.0f}"
                  f"{d['median_ha']:>11.2f}{d['median_shape']:>8.2f}")
        print("\n  observed area [ha] and ratio, by purity threshold:")
        print("  " + f"{'category':<12}" + "".join(f"{'>' + f'{t:g}':>16}"
                                                   for t in thresholds))
        for n in names:
            cells = "".join(f"{out[n]['observed_ha'][t]:>9,.0f}"
                            f"{out[n]['ratio'][t]:>7.2f}" for t in thresholds)
            print(f"  {n:<12}{cells}")
        print("\n  the ratio at >0 must exceed 1 — the PSF spreads past every boundary.")
    out["_window"] = win
    return out
