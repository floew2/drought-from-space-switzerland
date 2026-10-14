"""Reading GEE exports.

Two export generations are handled. Newer files carry `_median` and `_count`
per band; older ones only the median. The count matters: without it the
per-parcel coverage filter has nothing to work with and a parcel-date survives
on a single valid pixel.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .config import gee_dir

INDEX_COLS = ("parcel_id", "area_ha", "elev_m", "lnf_code")
INDICES = ("NDVI", "kNDVI", "NDMI", "NDRE")


def _is_wide(df: pd.DataFrame) -> bool:
    return any(c[0].isdigit() and ("NDVI" in c or "NDMI" in c)
               for c in df.columns)


def reshape_wide(df: pd.DataFrame) -> pd.DataFrame:
    """Wide GEE export -> tidy long frame.

    Band names come from the granule product ID, so a column looks like
        20260304T102921_20260304T103051_T32TMT_NDVI_median
        <sensing datetime>_<processing datetime>_<MGRS tile>_<index>_<stat>
    Only the first 8 characters (the sensing date), the index name and the
    statistic matter. The wide shape exists because ONE reduceRegions over a
    toBands() stack avoids the per-call overhead of one reduction per date —
    the difference between an affordable extraction and an unaffordable one.
    """
    id_cols = [c for c in INDEX_COLS if c in df.columns]
    val_cols = [c for c in df.columns
                if c[0].isdigit()
                and (c.endswith(("_NDVI", "_NDMI")) or "_NDVI_" in c
                     or "_NDMI_" in c)]
    if not val_cols:
        raise ValueError("no <date>_<tile>_<index> columns found")

    long = df.melt(id_vars=id_cols, value_vars=val_cols,
                   var_name="band", value_name="value").dropna(subset=["value"])

    has_stat = long["band"].str.endswith(("_median", "_count"))
    long["stat"] = np.where(has_stat,
                            long["band"].str.rsplit("_", n=1).str[-1], "median")
    stem = pd.Series(np.where(has_stat,
                              long["band"].str.rsplit("_", n=1).str[0],
                              long["band"]), index=long.index)
    long["index_name"] = stem.str.rsplit("_", n=1).str[-1]
    long["date"] = pd.to_datetime(stem.str.slice(0, 8), format="%Y%m%d")
    long["year"] = long["date"].dt.year
    long["doy"] = long["date"].dt.dayofyear

    keys = id_cols + ["year", "doy"]
    # Two MGRS tiles can see the same parcel on the same day. They observe the
    # SAME ground, so averaging is right (and a no-op when only one has data).
    med = (long[long.stat == "median"]
           .pivot_table(index=keys, columns="index_name", values="value",
                        aggfunc="mean").reset_index())
    med.columns.name = None

    if (long.stat == "count").any():
        cnt = (long[long.stat == "count"].groupby(keys)["value"].max()
               .rename("n_px").reset_index())
        med = med.merge(cnt, on=keys, how="left")
    else:
        med["n_px"] = 1.0
        med.attrs["no_pixel_count"] = True

    if "NDVI" in med:
        med["kNDVI"] = np.tanh(med["NDVI"] ** 2)   # a transform, not new data
    return med


def load_cuts(site: str, years=None, verbose: bool = True,
              labelled: bool = False) -> pd.DataFrame:
    """All per-parcel cut CSVs for a site (or only the labelled parcel-years)."""
    folder = gee_dir(site) / "labelled" if labelled else gee_dir(site)
    files = sorted(folder.glob("cuts_*.csv"))
    if years:
        files = [f for f in files
                 if any(str(y) in f.stem for y in years)]
    if not files:
        raise FileNotFoundError(
            f"no cuts_*.csv in {folder} — run `drought ee cuts"
            f"{' --labelled-only' if labelled else ''}`")

    frames, no_count = [], False
    from .provenance import verify
    for f in files:
        verify(f, "cuts")                 # refuses a stale export, naming why
        d = pd.read_csv(f)
        d.columns = [c.strip() for c in d.columns]
        if _is_wide(d):
            w = reshape_wide(d)
            no_count |= bool(w.attrs.get("no_pixel_count"))
            frames.append(w)
        else:
            ren = {c: c.replace("_median", "")
                   for c in d.columns if c.endswith("_median")}
            frames.append(d.rename(columns=ren)
                           .rename(columns={"NDVI_count": "n_px"}))
    df = pd.concat(frames, ignore_index=True)
    for c in ("year", "doy", "n_px", "area_ha", "elev_m", *INDICES):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["NDVI", "NDMI", "doy", "parcel_id"])

    # Safeguard, a no-op on the current export format: if a date appears twice
    # for one parcel, weight the index by pixel count but take n_px as the MAX.
    # Summing would double-count overlapping tiles, and one date with two clear
    # tiles would then set an inflated maximum that makes every ordinary date
    # look half-covered.
    def collapse(g):
        w = g["n_px"].to_numpy(float)
        out = {c: float(np.average(g[c], weights=w))
               for c in INDICES if c in g}
        out["n_px"] = float(w.max())
        return pd.Series(out)

    keys = ["parcel_id", "year", "doy"]
    meta_cols = [c for c in ("area_ha", "elev_m", "lnf_code") if c in df.columns]
    meta = df.groupby(keys)[meta_cols].first()
    df = (df.groupby(keys, as_index=False)
            .apply(collapse, include_groups=False).reset_index(drop=True))
    df = df.join(meta, on=keys)
    df.attrs["site"] = site
    df.attrs["no_pixel_count"] = no_count

    if verbose:
        print(f"{site}: {len(df):,} parcel-dates, "
              f"{df.parcel_id.nunique()} parcels, "
              f"years {sorted(df.year.unique())}")
        if no_count:
            print("  [!] this export carries no pixel count — the per-parcel "
                  "coverage filter cannot be applied. Re-extract to enable it.")
    return df


def load_trajectory(site: str | None = None, path=None) -> pd.DataFrame:
    """Site-mean index trajectories (one row per site/regime/year/date)."""
    if path is None:
        files = sorted(gee_dir(site).glob("*trajectory*.csv")) if site else []
        if not files:
            raise FileNotFoundError(
                f"no trajectory*.csv in {gee_dir(site) if site else '?'}")
        path = files[0]
    df = pd.read_csv(path)
    df.columns = [c.strip().lower() for c in df.columns]
    df = df.rename(columns={c.lower(): c for c in INDICES})
    for c in ("year", "doy", "n_px", *INDICES):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df
