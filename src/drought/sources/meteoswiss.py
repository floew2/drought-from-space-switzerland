"""MeteoSwiss Open Data: the SwissMetNet daily station series.

    Licence: CC BY 4.0, "Source: MeteoSwiss".
    Docs:    https://opendatadocs.meteoswiss.ch/a-data-groundbased/a1-automatic-weather-stations

`download` writes the station's daily files, as MeteoSwiss ships them, into
paths.station and records their provenance; `load_local` reads them back.
"""

from __future__ import annotations

import io
import logging
from pathlib import Path

import pandas as pd
import requests

__all__ = ["P_PRECIP", "P_TMAX", "P_TMEAN", "P_TMIN", "download", "load_local",
           "resolve_station"]

log = logging.getLogger(__name__)
TIMEOUT = 60

STAC = "https://data.geo.admin.ch/api/stac/v1/collections/ch.meteoschweiz.ogd-smn"
META_STATIONS = "https://data.geo.admin.ch/ch.meteoschweiz.ogd-smn/ogd-smn_meta_stations.csv"

# MeteoSwiss parameter short names (daily granularity)
P_PRECIP = "rre150d0"   # precipitation sum, 00:00-00:00 UTC+1 [mm]
P_TMEAN = "tre200d0"    # air temperature 2 m, daily mean [degC]
P_TMAX = "tre200dx"     # air temperature 2 m, daily maximum [degC]
P_TMIN = "tre200dn"     # air temperature 2 m, daily minimum [degC]


def _get(url: str) -> requests.Response:
    r = requests.get(url, timeout=TIMEOUT)
    r.raise_for_status()
    return r


def resolve_station(name: str, abbr: str | None = None) -> tuple[str, float, float]:
    """Return (abbr, latitude, altitude_m). Latitude is always resolved, because
    the PET calculation needs it even when the abbreviation is given explicitly."""
    meta = pd.read_csv(io.StringIO(_get(META_STATIONS).text), sep=";")
    if abbr:
        hit = meta[meta["station_abbr"].str.upper() == abbr.upper()]
    else:
        hit = meta[meta["station_name"].str.contains(name, case=False, na=False)]
    if hit.empty:
        raise LookupError(f"No station matching {abbr or name!r} in the MeteoSwiss station list.")
    if len(hit) > 1:
        log.info("Multiple matches, using the first:\n%s",
                 hit[["station_abbr", "station_name", "station_height_masl"]].to_string(index=False))
    row = hit.iloc[0]
    log.info("Station: %s (%s), %.0f m a.s.l., data since %s", row["station_name"],
             row["station_abbr"], row["station_height_masl"], row["station_data_since"])
    return (
        str(row["station_abbr"]).lower(),
        float(row["station_coordinates_wgs84_lat"]),
        float(row["station_height_masl"]),
    )


def daily_asset_urls(abbr: str) -> list[str]:
    """Discover the daily CSV assets via STAC instead of guessing filenames."""
    item = _get(f"{STAC}/items/{abbr}").json()
    urls = [
        a["href"]
        for name, a in item["assets"].items()
        if f"_{abbr}_d_" in name and name.endswith(".csv")
    ]
    if not urls:
        raise LookupError(f"No daily assets found for station {abbr!r}.")
    return sorted(urls)


def download(folder: str | Path, name: str, abbr: str | None = None) -> list[Path]:
    """Write a station's daily files and the station list into `folder`.

    Args:
        folder: Target directory (paths.station).
        name: Station name to look up, e.g. "St. Gallen".
        abbr: Station abbreviation, to skip the lookup.

    Returns:
        The daily files written, each with its provenance record.
    """
    from drought.core.provenance import register

    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "ogd-smn_meta_stations.csv").write_text(_get(META_STATIONS).text)
    abbr, _lat, _alt = resolve_station(name, abbr)
    out = []
    for url in daily_asset_urls(abbr):
        f = folder / url.rsplit("/", 1)[-1]
        f.write_text(_get(url).text)
        register(f, "station")
        out.append(f)
    log.info("  %d daily file(s) for %s -> %s", len(out), abbr, folder)
    return out


def load_local(folder, abbr: str | None = None) -> tuple[pd.DataFrame, str, float]:
    """Read daily files already on disk, plus the station latitude.

    Same layout as the MeteoSwiss Open Data download (one or more
    `ogd-smn_<abbr>_d_*.csv` and `ogd-smn_meta_stations.csv`), so a cached
    copy and the synthetic stand-in go through the same reader.

    Args:
        folder: Directory holding the files.
        abbr: Station abbreviation; inferred if exactly one station is present.

    Returns:
        (daily frame, abbr, latitude).

    Raises:
        FileNotFoundError: If no daily files are found.
    """
    folder = Path(folder)
    files = sorted(folder.glob(f"ogd-smn_{abbr or '*'}_d_*.csv"))
    if not files:
        raise FileNotFoundError(f"no ogd-smn_*_d_*.csv in {folder}")
    abbr = abbr or files[0].name.split("_")[1]
    meta = pd.read_csv(folder / "ogd-smn_meta_stations.csv", sep=";")
    row = meta[meta["station_abbr"].str.lower() == abbr.lower()].iloc[0]
    from drought.core.provenance import verify
    for f in files:
        verify(f, "station")
    return _read_daily(files), abbr.lower(), float(row["station_coordinates_wgs84_lat"])


def _read_daily(sources) -> pd.DataFrame:
    """Concatenate daily files (paths or buffers) into one typed frame."""
    frames = []
    for src in sources:
        url = getattr(src, "name", "buffer")
        df = pd.read_csv(src, sep=";")
        # MeteoSwiss ships inconsistent capitalisation between historical/recent
        df.columns = [c.strip().lower() for c in df.columns]
        ts_col = next(
            (c for c in ("reference_timestamp", "reference_ts", "date") if c in df.columns),
            None,
        )
        if ts_col is None:
            raise KeyError(f"No timestamp column in {url} (got {list(df.columns)[:6]})")
        df["date"] = pd.to_datetime(df[ts_col], format="%d.%m.%Y %H:%M", errors="coerce")
        if df["date"].isna().all():  # fall back if the format ever changes
            df["date"] = pd.to_datetime(df[ts_col], dayfirst=True, errors="coerce")
        frames.append(df)

    out = pd.concat(frames, ignore_index=True)
    out = out.dropna(subset=["date"]).drop_duplicates(subset="date", keep="last")
    for p in (P_PRECIP, P_TMEAN, P_TMAX, P_TMIN):
        if p in out.columns:
            out[p] = pd.to_numeric(out[p], errors="coerce")
    return out.set_index("date").sort_index()


