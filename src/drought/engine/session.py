"""Connection to Earth Engine and the Cloud Storage bucket.

Both use the credentials `earthengine authenticate` stores. Their scopes cover
Cloud Storage (devstorage.full_control), so the bucket needs no second login.
Storage is reached through its JSON API with an authorised session, which keeps
the dependency to earthengine-api alone.
"""

from __future__ import annotations

import functools
from pathlib import Path
from urllib.parse import quote

from ..config import CFG

__all__ = ["asset_id", "asset_properties", "bucket", "download", "ee", "exists",
           "init", "list_objects", "object_name", "project", "upload"]

_GCS = "https://storage.googleapis.com"


def project() -> str:
    """The configured Earth Engine project, or a message saying where to set it."""
    p = CFG.gee.project
    if not p:
        raise SystemExit("gee.project is not set. Put it in config.local.yaml:\n"
                         "  gee: {project: <your-ee-project>, bucket: <bucket>}")
    return str(p)


def bucket() -> str:
    """The configured Cloud Storage bucket."""
    b = CFG.gee.bucket
    if not b:
        raise SystemExit("gee.bucket is not set. Put it in config.local.yaml next to "
                         "gee.project. The bucket must be in a project with billing "
                         "enabled; the free tier covers this volume.")
    return str(b)


def object_name(*parts: str) -> str:
    """An object path under gee.bucket_prefix."""
    return "/".join([str(CFG.gee.bucket_prefix).strip("/"), *parts])


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


@functools.lru_cache(maxsize=1)
def _http():
    from google.auth.transport.requests import AuthorizedSession, Request

    creds = ee().data.get_persistent_credentials()
    if not creds.valid:
        creds.refresh(Request())
    return AuthorizedSession(creds)


def _check(r, what: str):
    if r.status_code >= 400:
        raise RuntimeError(f"Cloud Storage {what}: HTTP {r.status_code} {r.text[:300]}")
    return r


def upload(local: str | Path, name: str, content_type: str = "application/octet-stream") -> str:
    """Upload a file to the bucket. Returns its gs:// URI."""
    b = bucket()
    with open(local, "rb") as f:
        _check(_http().post(f"{_GCS}/upload/storage/v1/b/{b}/o",
                            params={"uploadType": "media", "name": name},
                            headers={"Content-Type": content_type}, data=f), "upload")
    return f"gs://{b}/{name}"


def download(name: str, dest: str | Path) -> Path:
    """Download one object to `dest`."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    r = _check(_http().get(f"{_GCS}/storage/v1/b/{bucket()}/o/{quote(name, safe='')}",
                           params={"alt": "media"}, stream=True), f"download {name}")
    with open(dest, "wb") as f:
        for chunk in r.iter_content(1 << 20):
            f.write(chunk)
    return dest


def list_objects(prefix: str) -> list[str]:
    """Object names under a prefix."""
    out, token = [], None
    while True:
        params = {"prefix": prefix, "fields": "items(name),nextPageToken"}
        if token:
            params["pageToken"] = token
        j = _check(_http().get(f"{_GCS}/storage/v1/b/{bucket()}/o", params=params),
                   "list").json()
        out += [i["name"] for i in j.get("items", [])]
        token = j.get("nextPageToken")
        if not token:
            return out


def exists(name: str) -> bool:
    """Is there an object of this name?"""
    return name in list_objects(name)


def asset_properties(name: str) -> dict | None:
    """Properties of an asset, or None if it does not exist."""
    try:
        return ee().data.getAsset(asset_id(name)).get("properties", {})
    except Exception as e:  # noqa: BLE001
        if "not found" in str(e).lower() or "does not exist" in str(e).lower():
            return None
        raise
