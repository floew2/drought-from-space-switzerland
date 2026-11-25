"""Provenance: every input and output carries a record of what produced it.

The failure this exists for: files repeatedly outlived the
method that made them. A CSV proposed for publication had the right schema
but the wrong vintage. Its management categories were ones the pipeline no
longer emits, its 2026 series ended before the result it was meant to support,
and the mask beside it was the old 20 m one. Nothing in the code could tell.
Being careful did not fix that; a structural check does:

* Each file gets a sidecar `<file>.provenance.json`: the settings that shape
  it, a fingerprint of the mask it depends on, a hash of the file itself, a
  summary of its content, the package version and a timestamp.
* `register` writes it when a file is produced, after checking that its
  content can be right under the current configuration.
* `verify` runs before a file is used. It refuses a missing record, a file
  changed since it was recorded, or any setting or mask that differs from the
  current ones, and it names the field that differs.

`config.yaml` `provenance.require: false` turns refusals into warnings.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import warnings
from pathlib import Path

from drought.core.config import CFG, ROOT, active_profile, mask_path

__all__ = ["KINDS", "ProvenanceError", "is_current", "mask_fingerprint", "record",
           "register", "settings", "settings_hash", "sidecar_path", "verify", "write_sidecar"]

# Settings that shape each kind of file. A change to any of them makes an older
# file of that kind stale. Paths and plotting options are deliberately absent.
_MASK = ("mask.res", "mask.subpixel", "mask.psf_model", "mask.psf_mtf_nyquist",
         "mask.sliver_floor_m2", "mask.min_overlap", "mask.min_years_same",
         "landuse.mown", "landuse.grazed", "landuse.intensive", "landuse.extensive")
_YEARS = ("project.current_year", "project.baseline_years")
# What decides which stored pixels are read, and how reflectance becomes an index.
_SOURCE = ("swisseo.collection", "swisseo.scale", "swisseo.offset", "swisseo.usable")
SETTING_KEYS = {
    "mask": _MASK,
    "condition": (*_MASK, *_YEARS, *_SOURCE, "mask.purity_min", "analysis.regimes",
                  "condition"),
    "sensitivity": (*_MASK, *_YEARS, *_SOURCE, "mask.purity_min", "analysis.regimes",
                    "sensitivity", "windows", "parcels.year"),
    "cuts": (*_MASK, *_YEARS, *_SOURCE, "mask.purity_min", "cuts", "parcels.year"),
    "anomaly": (*_MASK, *_YEARS, *_SOURCE, "mask.purity_min", "analysis.regimes",
                "anomaly", "windows", "parcels.year"),
    "vhi": (*_MASK, "project.current_year", "mask.purity_min", "analysis.regimes",
            "swisseo.vhi", "windows", "parcels.year"),
    "station": _YEARS,
    "results": (*_MASK, *_YEARS, *_SOURCE, "mask.purity_min", "analysis.regimes",
                "analysis.merge_regimes"),
}
MASK_DEPENDENT = {"condition", "cuts", "sensitivity", "anomaly", "vhi", "results"}
KINDS = tuple(SETTING_KEYS)


class ProvenanceError(RuntimeError):
    """A file's provenance does not match the current configuration.

    Attributes:
        path: The file.
        field: The first field that differs.
        recorded: Its value in the file's record.
        current: Its current value.
    """

    def __init__(self, path, field: str, recorded, current, hint: str = ""):
        self.path, self.field, self.recorded, self.current = path, field, recorded, current
        msg = (f"{Path(path).name}: {field} differs — recorded {recorded!r}, "
               f"current {current!r}.")
        super().__init__(msg + (f" {hint}" if hint else ""))


# --------------------------------------------------------------------------- #
# What a record contains
# --------------------------------------------------------------------------- #
def settings(kind: str = "results") -> dict:
    """The current values of the settings that shape a file of `kind`."""
    out = {}
    for key in SETTING_KEYS[kind]:
        v = CFG.get(key)
        out[key] = v.to_dict() if hasattr(v, "to_dict") else v
    return out


def settings_hash(values: dict | None = None) -> str:
    """Short, stable hash of a settings dict (default: the current settings)."""
    blob = json.dumps(values if values is not None else settings(), sort_keys=True)
    return hashlib.sha256(blob.encode()).hexdigest()[:12]


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def mask_fingerprint(path: str | Path | None = None) -> dict | None:
    """Identify a mask raster: content hash, resolution, bands, dtype, mtime.

    Args:
        path: The GeoTIFF (default: the configured mask).

    Returns:
        dict, or None if the mask does not exist.
    """
    import rasterio

    path = Path(path) if path else mask_path()
    if not path.exists():
        return None
    with rasterio.open(path) as ds:
        res = float(ds.res[0])
        bands = list(ds.descriptions)
        dtype = ds.dtypes[0]
    return {"path": _rel(path), "sha256": _sha(path), "res_m": res, "bands": bands,
            "dtype": dtype,
            "mtime": dt.datetime.fromtimestamp(path.stat().st_mtime).isoformat(timespec="seconds")}


def _rel(path: Path) -> str:
    try:
        return str(Path(path).resolve().relative_to(ROOT))
    except ValueError:
        return str(path)


def _version() -> str:
    from importlib.metadata import PackageNotFoundError, version
    try:
        return version("drought-from-space")
    except PackageNotFoundError:
        return "unknown"


def record(kind: str, path: str | Path | None = None, content: dict | None = None,
           depends_on_mask: bool | None = None) -> dict:
    """A provenance record for a file of `kind`, under the current configuration.

    Args:
        kind: What the file is, e.g. "trajectory", "cuts", "mask", "results".
        path: The file, hashed into the record if given.
        content: Summary of what is in the file (regimes, years, date range…).
        depends_on_mask: Include the current mask's fingerprint (default: by kind).

    Returns:
        The record as a dict.
    """
    if depends_on_mask is None:
        depends_on_mask = kind in MASK_DEPENDENT
    s = settings(kind)
    return {
        "kind": kind,
        "file": _rel(Path(path)) if path else None,
        "file_sha256": _sha(Path(path)) if path else None,
        "content": content or {},
        "settings": s,
        "settings_hash": settings_hash(s),
        "mask": mask_fingerprint() if depends_on_mask else None,
        "profile": active_profile(),
        "package_version": _version(),
        "created": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
    }


def sidecar_path(path: str | Path) -> Path:
    """`data.csv` -> `data.csv.provenance.json`."""
    path = Path(path)
    return path.with_name(path.name + ".provenance.json")


def write_sidecar(path: str | Path, rec: dict) -> Path:
    """Write a record beside its file. Returns the sidecar path."""
    out = sidecar_path(path)
    out.write_text(json.dumps(rec, indent=1, default=str))
    return out


# --------------------------------------------------------------------------- #
# Content summaries, and what they must satisfy
# --------------------------------------------------------------------------- #
def summarise(kind: str, path: str | Path) -> dict:
    """What a data file contains, read from the file itself."""
    import pandas as pd

    if kind == "station":
        st = pd.read_csv(path, sep=";")
        ts = pd.to_datetime(st.iloc[:, 1], format="%d.%m.%Y %H:%M", errors="coerce")
        return {"rows": len(st), "first_day": str(ts.min().date()),
                "last_day": str(ts.max().date())}
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    if kind == "condition":
        cur = int(CFG.project.current_year)
        last = df.loc[df["year"] == cur, "doy"].max() if "year" in df else None
        return {"rows": len(df), "sites": sorted(df["site"].unique().tolist()),
                "regimes": sorted(df["regime"].unique().tolist()),
                "years": sorted(int(y) for y in df["year"].unique()),
                "current_year_last_doy": None if last is None or last != last else int(last)}
    if kind == "cuts":
        dates = sorted({c[:8] for c in df.columns if c[:1].isdigit()})
        return {"parcels": len(df), "dates": len(dates),
                "first_date": dates[0] if dates else None,
                "last_date": dates[-1] if dates else None,
                "has_geometry": ".geo" in df.columns}
    if kind == "sensitivity":
        return {"rows": len(df), "years": sorted(int(y) for y in df["year"].unique()),
                "classes": int(df["pclass"].nunique())}
    return {"rows": len(df)}


def _check_content(kind: str, path, content: dict) -> None:
    """Content that cannot be right under the current configuration."""
    if kind == "condition":
        expected = sorted(CFG.analysis.regimes)
        if content.get("regimes") != expected:
            raise ProvenanceError(path, "content.regimes", content.get("regimes"), expected,
                                  "The series carry management categories the pipeline "
                                  "does not use: they predate the current method.")
        years = sorted([*CFG.project.baseline_years, CFG.project.current_year])
        missing = [y for y in years if y not in content.get("years", [])]
        if missing:
            raise ProvenanceError(path, "content.years", content.get("years"), years,
                                  f"Missing years {missing}.")
    if kind == "station":
        return
    if kind == "cuts" and content.get("has_geometry"):
        raise ProvenanceError(path, "content.has_geometry", True, False,
                              "A cut series with a geometry column must not enter the "
                              "repository; strip it before registering.")


# --------------------------------------------------------------------------- #
# Registering, and verifying before use
# --------------------------------------------------------------------------- #
def register(path: str | Path, kind: str) -> Path:
    """Record a file as it is produced, after checking it can be current.

    Args:
        path: The file.
        kind: One of KINDS.

    Returns:
        The sidecar path.

    Raises:
        ProvenanceError: If its content cannot be right under the current
            configuration (for example the wrong management categories).
    """
    path = Path(path)
    content = summarise(kind, path)
    _check_content(kind, path, content)
    return write_sidecar(path, record(kind, path, content))


def verify(path: str | Path, kind: str) -> dict:
    """Check a file against the current configuration before using it.

    Args:
        path: The file.
        kind: Expected kind of record.

    Returns:
        The file's record.

    Raises:
        ProvenanceError: Naming the first field that differs (unless
            `provenance.require` is false, in which case it warns).
    """
    try:
        return _verify(Path(path), kind)
    except ProvenanceError as exc:
        if getattr(CFG, "provenance", None) is not None and not CFG.provenance.require:
            warnings.warn(f"provenance not verified: {exc}", stacklevel=2)
            return {}
        raise


def _verify(path: Path, kind: str) -> dict:
    side = sidecar_path(path)
    if not side.exists():
        raise ProvenanceError(path, "provenance record", None, "required",
                              f"No {side.name}: the file was not made by this pipeline, "
                              "or its record was deleted. Re-run the stage that makes it.")
    rec = json.loads(side.read_text())
    if rec.get("kind") != kind:
        raise ProvenanceError(path, "kind", rec.get("kind"), kind)
    if rec.get("file_sha256") and rec["file_sha256"] != _sha(path):
        raise ProvenanceError(path, "file_sha256", rec["file_sha256"], _sha(path),
                              "The file changed after it was recorded.")
    now = settings(kind)
    for key in SETTING_KEYS[kind]:
        if rec["settings"].get(key) != now.get(key):
            raise ProvenanceError(path, key, rec["settings"].get(key), now.get(key),
                                  "The file was produced under different settings.")
    if kind == "mask":
        from drought.core.config import CFG as _C
        fp = mask_fingerprint(path)
        if fp is None:
            raise ProvenanceError(path, "mask", None, None, "The mask file is missing.")
        if fp["res_m"] != float(_C.mask.res) or fp["dtype"] != "float32":
            raise ProvenanceError(path, "mask resolution / dtype",
                                  f"{fp['res_m']:g} m {fp['dtype']}",
                                  f"{float(_C.mask.res):g} m float32",
                                  "This is not a mask the current pipeline builds.")
    if rec.get("mask") is not None:
        cur = mask_fingerprint()
        if cur is None:
            raise ProvenanceError(path, "mask", rec["mask"]["sha256"], None,
                                  "The mask it depends on is not built.")
        for f in ("res_m", "bands", "dtype", "sha256"):
            if rec["mask"].get(f) != cur.get(f):
                raise ProvenanceError(path, f"mask.{f}", rec["mask"].get(f), cur.get(f),
                                      "It was produced with a different mask.")
    _check_content(kind, path, rec.get("content", {}))
    return rec


def is_current(path: str | Path, kind: str) -> bool:
    """Whether a file exists and passes `verify`, without raising or warning."""
    path = Path(path)
    if not path.exists():
        return False
    try:
        _verify(path, kind)
    except ProvenanceError:
        return False
    return True
