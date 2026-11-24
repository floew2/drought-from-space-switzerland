"""Trajectory, purity sensitivity and cut series from swissEO, in one pass.

Each mosaic is read once per site, and all three products are taken from it:

trajectory   per management system, the mean of the four indices over the
             cells whose purity reaches mask.purity_min, and their count; kept
             if the count reaches gee.trajectory.min_pixels. Same columns as
             the Earth Engine export.
sensitivity  per year, the per-cell median NDMI over the site's window
             (gee.sensitivity.windows), grouped by purity class, over the cells
             of parcels of that system (as in Earth Engine: inside the parcel
             outline, purity > 0, parcels with a cell at gee.strict_purity).
cuts         per parcel, the median and count of NDVI and NDMI over its pure
             mown cells after the edge buffer, one column per mosaic. The
             random sample (gee.cuts.sample_n, seed) and the hand-labelled
             parcels are taken in the same pass.

A mosaic whose share of not-clear grassland cells exceeds swisseo.max_site_cloud
is left out of the sensitivity composite and the cut series, as Earth Engine
leaves out granules above gee.*.max_scene_cloud; the trajectory keeps it if
enough cells are clear, as in Earth Engine.

Parcel ids. Earth Engine names a feature of an uploaded shapefile by its
position in the file, as 20 hexadecimal digits; the hand labels carry those ids.
Here a parcel's id is built the same way from the same file
(paths.parcels, gee.parcels_year). `drought swisseo check-ids` confirms that
both sides name the same fields before labels are scored on swissEO series.

Results per mosaic are cached under <paths.swisseo>/_cache/<CONFIG_HASH>/, so a
run that stops resumes, and a second run costs nothing.
"""

from __future__ import annotations

import datetime as dt
import hashlib
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from ..config import CFG, ROOT, series_dir, series_root
from . import grid, scene, stac

__all__ = ["Plan", "parcels_table", "plan", "process", "run", "years"]

TRAJ_COLS = ["site", "regime", "year", "date", "doy", "NDVI", "kNDVI", "NDMI", "NDRE", "n_px"]
SENS_COLS = ["year", "pclass", "purity_lo", "mean", "n_px"]


def years() -> list[int]:
    return [*CFG.project.baseline_years, CFG.project.current_year]


def config_hash() -> str:
    from ..gee import _remember, expected_hash
    h = expected_hash(with_mask=True)
    _remember(ROOT / "output" / "gee", h)     # so a mismatch can name what changed
    return h


# --------------------------------------------------------------------------- #
# Parcels
# --------------------------------------------------------------------------- #
def parcels_table(path: Path | None = None):
    """The cadastre the per-parcel work runs on, with Earth Engine's ids.

    Returns:
        GeoDataFrame in LV95: parcel_id (20 hex digits of the row position),
        lnf_code, area_ha (planar, before the edge buffer), geometry.
    """
    import geopandas as gpd

    from ..mask.parcels import _year_files
    from ..mask.settings import CODE_FIELD_CANDIDATES

    path = path or _year_files()[int(CFG.gee.parcels_year)]
    g = gpd.read_file(path)
    if g.crs is None:
        raise ValueError(f"{path}: no CRS")
    g = g.to_crs(2056)
    code = next((c for c in CODE_FIELD_CANDIDATES if c in g.columns), None)
    if code is None:
        raise ValueError(f"{path}: none of {CODE_FIELD_CANDIDATES} in {list(g.columns)}")
    out = gpd.GeoDataFrame({
        "parcel_id": [f"{i:020x}" for i in range(len(g))],
        "lnf_code": pd.to_numeric(g[code], errors="coerce"),
        "area_ha": g.geometry.area / 1e4,
    }, geometry=g.geometry.values, crs=2056)
    return out[out.geometry.notna() & ~out.geometry.is_empty]


def _elevation(geoms) -> np.ndarray:
    """DEM value at each parcel's representative point, or NaN without a DEM."""
    dem = getattr(CFG.swisseo, "dem", None)
    if not dem:
        return np.full(len(geoms), np.nan)
    import rasterio
    pts = [(p.x, p.y) for p in geoms.representative_point()]
    with rasterio.open(ROOT / dem) as ds:
        v = np.array([x[0] for x in ds.sample(pts)], dtype=float)
        if ds.nodata is not None:
            v[v == ds.nodata] = np.nan
    return v


def _labelled(site: str) -> dict[int, list[str]]:
    from ..labels import load_labels
    out: dict[int, list[str]] = {}
    try:
        keys = load_labels(site)
    except FileNotFoundError:
        return out
    for key in keys:
        pid, year = key.rsplit("|", 1)
        out.setdefault(int(year), []).append(pid)
    return {y: sorted(set(v)) for y, v in out.items()}


# --------------------------------------------------------------------------- #
# What is measured where, per site
# --------------------------------------------------------------------------- #
@dataclass
class Plan:
    site: str
    win: grid.Window
    regimes: dict[str, np.ndarray]          # regime -> pure cells (bool)
    grass: np.ndarray                       # any group purity > 0: the cloud test
    sens_cells: dict[str, np.ndarray]       # regime -> flat indices of sensitivity cells
    sens_class: dict[str, np.ndarray]       # regime -> purity class of those cells
    parcel_label: np.ndarray                # 0 none, k = row k-1 of `parcels`
    parcels: pd.DataFrame                   # parcel_id, area_ha, elev_m, lnf_code
    sample_ids: list[str]
    labelled: dict[int, list[str]] = field(default_factory=dict)

    @property
    def parcel_digest(self) -> str:
        return hashlib.sha256("|".join(self.parcels.parcel_id).encode()).hexdigest()[:12]


def _rasterize(geoms, values, win: grid.Window) -> np.ndarray:
    """Cell centres inside each polygon, as Earth Engine reduces a polygon."""
    from rasterio.features import rasterize
    shapes = [(g, int(v)) for g, v in zip(geoms, values) if g is not None and not g.is_empty]
    if not shapes:
        return np.zeros(win.shape, np.int32)
    return rasterize(shapes, out_shape=win.shape, transform=win.transform,
                     fill=0, all_touched=False, dtype="int32")


def plan(site: str, parcels=None, mask_file=None) -> Plan:
    """Cells, classes and parcels for one site, from the mask and the cadastre."""
    from shapely.geometry import box

    win = grid.site_window(site, mask_file)
    m = grid.read_mask(win, mask_file)
    pmin = float(CFG.mask.purity_min)
    regimes = {r: m[r] >= pmin for r in CFG.analysis.regimes}
    grass = np.zeros(win.shape, bool)
    for a in m.values():
        grass |= a > 0

    p = parcels_table() if parcels is None else parcels
    p = p[p.intersects(box(*win.bounds))].reset_index(drop=True)
    lu = CFG.landuse

    # Sensitivity: cells inside parcels of the system, purity > 0, and parcels
    # holding a cell at strict_purity (0.0 = every parcel).
    n = int(CFG.gee.sensitivity.n_classes)
    strict = float(CFG.gee.strict_purity)
    sens_cells, sens_class = {}, {}
    for r in CFG.analysis.regimes:
        sub = p[p.lnf_code.isin(list(getattr(lu, r)))]
        lab = _rasterize(sub.geometry, np.arange(1, len(sub) + 1), win)
        if strict > 0:
            hit = np.unique(lab[(lab > 0) & (m[r] >= strict)])
            lab = np.where(np.isin(lab, hit), lab, 0)
        cells = (lab > 0) & (m[r] > 0)
        idx = np.flatnonzero(cells)
        sens_cells[r] = idx
        sens_class[r] = np.minimum(np.floor(m[r].ravel()[idx] * n), n - 1).astype(np.int16)

    # Cuts: mown parcels of at least min_area_ha, shrunk by the edge buffer,
    # measured on their pure mown cells.
    c = CFG.gee.cuts
    mown = p[p.lnf_code.isin(list(lu.mown)) & (p.area_ha >= float(c.min_area_ha))]
    mown = mown.sort_values("parcel_id").reset_index(drop=True)
    rng = np.random.default_rng(int(c.seed))
    k = min(int(c.sample_n), len(mown))
    sample_ids = sorted(mown.parcel_id.iloc[rng.choice(len(mown), k, replace=False)]) if k else []
    labelled = _labelled(site)
    wanted = set(sample_ids) | {i for v in labelled.values() for i in v}
    chosen = p[p.parcel_id.isin(wanted)].sort_values("parcel_id").reset_index(drop=True)
    missing = wanted - set(chosen.parcel_id)
    if missing:
        print(f"  [!] {site}: {len(missing)} labelled parcel id(s) not in the cadastre "
              f"inside the site window, e.g. {sorted(missing)[:3]} — "
              "run `drought swisseo check-ids`")
    shrunk = chosen.geometry.buffer(float(c.edge_buffer_m))
    lab = _rasterize(shrunk, np.arange(1, len(chosen) + 1), win)
    lab[~(m["mown"] >= pmin)] = 0
    table = pd.DataFrame({"parcel_id": chosen.parcel_id, "area_ha": chosen.area_ha,
                          "elev_m": _elevation(chosen.geometry),
                          "lnf_code": chosen.lnf_code.astype("Int64")})
    if table.elev_m.isna().all():
        # Without a DEM there is no elevation. Leave the column out rather than
        # empty: a key column of NaN would drop every row when reshaped.
        table = table.drop(columns="elev_m")
    return Plan(site, win, regimes, grass, sens_cells, sens_class, lab, table,
                sample_ids, labelled)


# --------------------------------------------------------------------------- #
# One mosaic
# --------------------------------------------------------------------------- #
def _in_window(doy: int, w) -> bool:
    return int(w[0]) <= doy <= int(w[1])


def process(item: stac.Item, pl: Plan) -> dict:
    """Everything the three products need from one mosaic, small enough to cache."""
    ok = scene.clear_mask(item, pl.win)
    n_grass = int(pl.grass.sum())
    cloud_pct = 100.0 * (1 - (ok & pl.grass).sum() / n_grass) if n_grass else 100.0
    out = {"cloud_pct": cloud_pct, "traj": {}, "sens": {}, "cuts": None}
    usable = ok & pl.grass
    if not usable.any():
        return out
    ix = scene.indices(item, pl.win, ok)

    for r, cells in pl.regimes.items():
        sel = cells & ~np.isnan(ix["NDVI"])
        n = int(sel.sum())
        if n:
            out["traj"][r] = [float(np.nanmean(ix[i][cells])) for i in scene.INDICES] + [n]

    if cloud_pct <= float(CFG.swisseo.max_site_cloud):
        w = getattr(CFG.gee.sensitivity.windows, pl.site, None)
        if w is not None and _in_window(item.doy, w):
            flat = ix[CFG.gee.sensitivity.index].ravel()
            out["sens"] = {r: flat[c].astype("float16") for r, c in pl.sens_cells.items()}
        lab = pl.parcel_label.ravel()
        has = lab > 0
        if has.any():
            d = pd.DataFrame({"k": lab[has], "NDVI": ix["NDVI"].ravel()[has],
                              "NDMI": ix["NDMI"].ravel()[has]})
            g = d.groupby("k").agg(NDVI_median=("NDVI", "median"), NDVI_count=("NDVI", "count"),
                                   NDMI_median=("NDMI", "median"), NDMI_count=("NDMI", "count"))
            out["cuts"] = g.reindex(np.arange(1, len(pl.parcels) + 1))
    return out


def _cache_file(root: Path, item: stac.Item, pl: Plan) -> Path:
    return root / pl.site / f"{item.id}_{pl.parcel_digest}.npz"


def _save(path: Path, r: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    arrs = {"cloud_pct": np.array(r["cloud_pct"])}
    for k, v in r["traj"].items():
        arrs[f"traj_{k}"] = np.array(v, dtype=float)
    for k, v in r["sens"].items():
        arrs[f"sens_{k}"] = v
    if r["cuts"] is not None:
        arrs["cuts"] = r["cuts"].to_numpy(dtype=float)
    tmp = path.with_suffix(".tmp.npz")
    np.savez_compressed(tmp, **arrs)
    tmp.replace(path)


def _load(path: Path) -> dict:
    z = np.load(path)
    out = {"cloud_pct": float(z["cloud_pct"]), "traj": {}, "sens": {}, "cuts": None}
    for k in z.files:
        if k.startswith("traj_"):
            out["traj"][k[5:]] = z[k].tolist()
        elif k.startswith("sens_"):
            out["sens"][k[5:]] = z[k]
        elif k == "cuts":
            out["cuts"] = pd.DataFrame(z[k], columns=["NDVI_median", "NDVI_count",
                                                      "NDMI_median", "NDMI_count"])
    return out


# --------------------------------------------------------------------------- #
# The run
# --------------------------------------------------------------------------- #
def _season() -> tuple[int, int]:
    ws = [CFG.gee.trajectory.doy, CFG.gee.cuts.doy,
          *CFG.gee.sensitivity.windows.to_dict().values()]
    return min(int(w[0]) for w in ws), max(int(w[1]) for w in ws)


def items_for(site: str, year: int, session=None) -> list[stac.Item]:
    d0, d1 = _season()
    start = dt.date(year, 1, 1) + dt.timedelta(days=d0 - 1)
    end = dt.date(year, 1, 1) + dt.timedelta(days=d1 - 1)
    return stac.search(getattr(CFG.sites, site).bbox_wgs84, start, end, session=session)


def _tag(item: stac.Item) -> str:
    return f"{item.datetime:%Y%m%dT%H%M%S}_SWISSEO"


def _write(df: pd.DataFrame, path: Path, kind: str) -> Path:
    from ..provenance import register
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    register(path, kind)
    rel = path.relative_to(ROOT) if path.is_relative_to(ROOT) else path
    print(f"  -> {rel}  (registered)")
    return path


def run(sites=None, only_years=None, workers: int = 4, dry_run: bool = False,
        session=None, parcels=None) -> list[Path]:
    """Read swissEO for every site and year; write and register the three products.

    Args:
        sites: Default: the union of gee.trajectory.sites and gee.cuts.sites.
        only_years: Restrict to these years. The trajectory and the sensitivity
            are written only when every year was read (their content check
            requires all years); the cut series are written per year.
        workers: Mosaics read in parallel.
        dry_run: List how many mosaics each site and year has, read nothing.
        session, parcels: For tests.

    Returns:
        The files written.
    """
    from ..engine import cuts as ee_cuts
    from ..engine import sensitivity as ee_sens
    from ..engine import trajectory as ee_traj

    h = config_hash()
    sites = sites or list(dict.fromkeys([*CFG.gee.trajectory.sites, *CFG.gee.cuts.sites]))
    yrs = sorted(only_years) if only_years else years()
    complete = set(yrs) == set(years())
    cache = series_root() / "_cache" / h
    traj_rows, written = [], []
    p_all = parcels

    for site in sites:
        if dry_run:
            for y in yrs:
                print(f"  {site} {y}: {len(items_for(site, y, session))} mosaics")
            continue
        if p_all is None:
            p_all = parcels_table()
        pl = plan(site, p_all)
        print(f"  {site}: window {pl.win.width} x {pl.win.height} cells; "
              + ", ".join(f"{r} {int(c.sum()):,} pure" for r, c in pl.regimes.items())
              + f"; {len(pl.sample_ids)} sampled + labelled parcels -> {len(pl.parcels)}")
        sens = {r: [] for r in pl.regimes}
        for year in yrs:
            items = items_for(site, year, session)

            def one(it, pl=pl):
                f = _cache_file(cache, it, pl)
                if f.exists():
                    return it, _load(f)
                with scene.gdal_env():
                    r = process(it, pl)
                _save(f, r)
                return it, r

            with ThreadPoolExecutor(max_workers=max(1, workers)) as ex:
                res = list(ex.map(one, items))
            print(f"  {site} {year}: {len(items)} mosaics, "
                  f"{sum(1 for _, r in res if r['traj'])} with clear grassland")

            # trajectory
            d0, d1 = CFG.gee.trajectory.doy
            if site in list(CFG.gee.trajectory.sites):
                for it, r in res:
                    if not (int(d0) <= it.doy <= int(d1)):
                        continue
                    for reg, v in r["traj"].items():
                        if v[-1] >= int(CFG.gee.trajectory.min_pixels) and \
                                not any(np.isnan(v[:4])):
                            traj_rows.append([site, reg, year, f"{it.datetime:%Y-%m-%d}",
                                              it.doy, *v[:4], int(v[4])])

            # sensitivity: per-cell median of the year, then class means
            for reg in pl.regimes:
                stack = [r["sens"][reg].astype("float32") for _, r in res if reg in r["sens"]]
                if not stack:
                    continue
                with np.errstate(all="ignore"):
                    import warnings
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", RuntimeWarning)
                        comp = np.nanmedian(np.vstack(stack), axis=0)
                d = pd.DataFrame({"pclass": pl.sens_class[reg], "v": comp}).dropna()
                g = d.groupby("pclass")["v"].agg(["mean", "count"]).reset_index()
                n = int(CFG.gee.sensitivity.n_classes)
                sens[reg].append(pd.DataFrame({
                    "year": year, "pclass": g.pclass, "purity_lo": g.pclass / n,
                    "mean": g["mean"], "n_px": g["count"]}))

            # cuts
            if site in list(CFG.gee.cuts.sites):
                c0, c1 = CFG.gee.cuts.doy
                cols = {}
                for it, r in res:
                    if r["cuts"] is None or not (int(c0) <= it.doy <= int(c1)):
                        continue
                    t = _tag(it)
                    for col in r["cuts"].columns:
                        name, stat = col.split("_")
                        cols[f"{t}_{name}_{stat}"] = r["cuts"][col].to_numpy()
                wide = pd.concat([pl.parcels.reset_index(drop=True),
                                  pd.DataFrame(cols)], axis=1)
                sample = wide[wide.parcel_id.isin(pl.sample_ids)]
                written.append(_write(sample, series_dir(site) /
                                      f"{ee_cuts.prefix(site, year, h)}.csv", "cuts"))
                ids = pl.labelled.get(year)
                if ids:
                    lab = wide[wide.parcel_id.isin(ids)]
                    written.append(_write(lab, series_dir(site) / "labelled" /
                                          f"cuts_{site}_{year}_labelled_cfg{h}.csv", "cuts"))

        if complete and site in CFG.gee.sensitivity.windows.to_dict():
            for reg, frames in sens.items():
                if frames:
                    df = pd.concat(frames, ignore_index=True)[SENS_COLS]
                    written.append(_write(df, series_root() /
                                          f"{ee_sens.prefix(site, reg, h)}.csv", "sensitivity"))

    if dry_run:
        return []
    if complete and traj_rows:
        df = pd.DataFrame(traj_rows, columns=TRAJ_COLS)
        written.insert(0, _write(df, series_root() / f"{ee_traj.prefix(h)}.csv", "trajectory"))
    elif not complete:
        print("  trajectory and sensitivity not written: they need every year "
              f"({years()}); the per-date results are cached for the full run")
    return written
