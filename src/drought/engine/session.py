"""Connection to Earth Engine, and downloads from Google Drive.

Exports go to Drive (batch tasks, no interactive limits, no billing), and are
fetched through the Drive API with the credentials `earthengine authenticate`
stores. Their default scopes include Drive, so there is no second login.
"""

from __future__ import annotations

import functools
from pathlib import Path

from ..config import CFG

__all__ = ["asset_id", "asset_info", "drive_download", "drive_find", "ee", "init",
           "project"]

_DRIVE = "https://www.googleapis.com/drive/v3/files"


def project() -> str:
    """The configured Earth Engine project, or a message saying where to set it."""
    p = CFG.gee.project
    if not p:
        raise SystemExit("gee.project is not set. Put it in config.local.yaml:\n"
                         "  gee: {project: <your-ee-project>}")
    return str(p)


def asset_id(name: str) -> str:
    """Full asset path in the configured project."""
    return f"projects/{project()}/assets/{name}"


@functools.lru_cache(maxsize=1)
def ee():
    """The initialised earthengine-api module."""
    try:
        import ee as _ee
    except ImportError as e:
        raise SystemExit("earthengine-api is not installed: uv sync --extra ee") from e
    try:
        _ee.Initialize(project=project())
    except Exception as e:  # noqa: BLE001 — ee raises several types here
        raise SystemExit(f"Earth Engine did not initialise ({e}).\n"
                         "Run `earthengine authenticate` once, then try again.") from e
    return _ee


def init():
    """Initialise Earth Engine; return the module."""
    return ee()


def asset_info(name: str) -> dict | None:
    """An asset's metadata (properties, bands, ...), or None if it does not exist."""
    try:
        return ee().data.getAsset(asset_id(name))
    except Exception as e:  # noqa: BLE001
        if "not found" in str(e).lower() or "does not exist" in str(e).lower():
            return None
        raise


# --------------------------------------------------------------------------- #
# Google Drive
# --------------------------------------------------------------------------- #
@functools.lru_cache(maxsize=2)
def _http(quota_project: bool):
    from google.auth.transport.requests import AuthorizedSession, Request

    creds = ee().data.get_persistent_credentials()
    if not quota_project and getattr(creds, "quota_project_id", None):
        creds = creds.with_quota_project(None)
    if not creds.valid:
        creds.refresh(Request())
    return AuthorizedSession(creds)


def _get(url: str, **kw):
    """GET against the Drive API, with a message that says what to fix."""
    r = _http(False).get(url, **kw)
    if r.status_code == 403:
        r = _http(True).get(url, **kw)          # with your project as quota project
    if r.status_code == 403:
        raise SystemExit(
            "Google Drive refused the request (HTTP 403). Either the Earth Engine "
            "login lacks Drive access (run `earthengine authenticate --force` and "
            "allow Drive), or the Drive API is off in your project (Cloud Console → "
            "APIs & Services → enable 'Google Drive API'; free, no billing).\n"
            + r.text[:300])
    if r.status_code >= 400:
        raise RuntimeError(f"Google Drive: HTTP {r.status_code} {r.text[:300]}")
    return r


def drive_find(filename: str) -> list[dict]:
    """Files of exactly this name in Drive, newest first."""
    q = f"name = '{filename}' and trashed = false"
    j = _get(_DRIVE, params={"q": q, "orderBy": "createdTime desc",
                             "fields": "files(id,name,createdTime,size)"}).json()
    return j.get("files", [])


def drive_download(file_id: str, dest: str | Path) -> Path:
    """Download one Drive file to `dest`."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    r = _get(f"{_DRIVE}/{file_id}", params={"alt": "media"}, stream=True)
    with open(dest, "wb") as f:
        for chunk in r.iter_content(1 << 20):
            f.write(chunk)
    return dest
