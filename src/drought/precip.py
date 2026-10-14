"""
Cumulative precipitation, St. Gallen: current year vs. 1991-2020 climatology.

Single source for BOTH the current year and the baseline: MeteoSwiss Open Data
(SwissMetNet automatic weather stations, daily values).
    Licence: CC-BY. When redistributing, credit "Quelle: MeteoSchweiz".
    Docs:    https://opendatadocs.meteoswiss.ch/a-data-groundbased/a1-automatic-weather-stations

Key design decisions (vs. a naive monthly-cumsum plot):
  * daily resolution, so the current year is truncated at exactly the same
    day-of-year as the baseline -> no partial-month artefact;
  * the baseline is derived from the same station and same instrument network,
    not hard-coded and not taken from a different source;
  * a 10th-90th percentile envelope of the 30 individual baseline years is shown,
    so the reader can see whether the current year is outside the historical
    range, not merely "below average";
  * optional second panel: cumulative climatic water balance P - PET
    (Hargreaves), because precipitation alone understates a hot drought.

Requires: pandas, numpy, matplotlib, requests
"""

from __future__ import annotations

import io
import textwrap

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests

# --------------------------------------------------------------------------- #
# Configuration — from config.yaml (precip: and project: sections)
# --------------------------------------------------------------------------- #
from .config import CFG, figures_dir

STATION_NAME = CFG.precip.station_name
STATION_ABBR = None              # set e.g. "STG" to skip the lookup
BASELINE = tuple(CFG.precip.baseline)
CURRENT_YEAR = CFG.project.current_year
INCLUDE_WATER_BALANCE = CFG.precip.include_water_balance
TIMEOUT = 60

STAC = "https://data.geo.admin.ch/api/stac/v1/collections/ch.meteoschweiz.ogd-smn"
META_STATIONS = "https://data.geo.admin.ch/ch.meteoschweiz.ogd-smn/ogd-smn_meta_stations.csv"

# MeteoSwiss parameter short names (daily granularity)
P_PRECIP = "rre150d0"   # precipitation sum, 00:00-00:00 UTC+1 [mm]
P_TMEAN = "tre200d0"    # air temperature 2 m, daily mean [degC]
P_TMAX = "tre200dx"     # air temperature 2 m, daily maximum [degC]
P_TMIN = "tre200dn"     # air temperature 2 m, daily minimum [degC]


# --------------------------------------------------------------------------- #
# Data access
# --------------------------------------------------------------------------- #
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
        print("Multiple matches, using the first:")
        print(hit[["station_abbr", "station_name", "station_height_masl"]].to_string(index=False))
    row = hit.iloc[0]
    print(
        f"Station: {row['station_name']} ({row['station_abbr']}), "
        f"{row['station_height_masl']:.0f} m a.s.l., "
        f"data since {row['station_data_since']}"
    )
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


def load_daily(abbr: str) -> pd.DataFrame:
    """Download all daily files for a station into one tidy, typed frame."""
    return _read_daily([io.StringIO(_get(url).text) for url in daily_asset_urls(abbr)])


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
    from pathlib import Path

    folder = Path(folder)
    files = sorted(folder.glob(f"ogd-smn_{abbr or '*'}_d_*.csv"))
    if not files:
        raise FileNotFoundError(f"no ogd-smn_*_d_*.csv in {folder}")
    abbr = abbr or files[0].name.split("_")[1]
    meta = pd.read_csv(folder / "ogd-smn_meta_stations.csv", sep=";")
    row = meta[meta["station_abbr"].str.lower() == abbr.lower()].iloc[0]
    from .provenance import verify
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


# --------------------------------------------------------------------------- #
# Day-of-year alignment
# --------------------------------------------------------------------------- #
def doy_matrix(series: pd.Series, years: range) -> pd.DataFrame:
    """Wide frame: rows = calendar-aligned DOY 1..365, cols = year.

    29 February is dropped and later DOYs shifted, so a given row means the same
    calendar date in leap and non-leap years alike.
    """
    df = series.to_frame("v")
    df["year"] = df.index.year
    df["doy"] = df.index.dayofyear
    leap = df.index.is_leap_year
    after_feb29 = df["doy"] > 60
    df = df[~(leap & (df["doy"] == 60))]                       # drop 29 Feb
    df.loc[leap & after_feb29, "doy"] -= 1                     # realign the rest
    wide = df.pivot_table(index="doy", columns="year", values="v", aggfunc="first")
    return wide.reindex(index=range(1, 366), columns=list(years))


def cumulative(wide: pd.DataFrame) -> pd.DataFrame:
    """Cumulative sum per year; NaN days are treated as 0 but reported upstream."""
    return wide.fillna(0.0).cumsum()


def report_gaps(wide: pd.DataFrame, label: str, upto: int | None = None) -> None:
    sub = wide.iloc[:upto] if upto else wide
    gaps = sub.isna().sum()
    bad = gaps[gaps > 0]
    if not bad.empty:
        print(f"  [!] {label}: missing days per year -> {bad.to_dict()}")
    else:
        print(f"  {label}: no missing days")


# --------------------------------------------------------------------------- #
# Hargreaves reference evapotranspiration
# --------------------------------------------------------------------------- #
def hargreaves_pet(tmean, tmax, tmin, doy, lat_deg):
    """Reference ET0 [mm/day], Hargreaves-Samani (FAO-56 Ch. 4 fallback method)."""
    phi = np.radians(lat_deg)
    j = np.asarray(doy, dtype=float)
    dr = 1 + 0.033 * np.cos(2 * np.pi * j / 365.0)
    delta = 0.409 * np.sin(2 * np.pi * j / 365.0 - 1.39)
    ws = np.arccos(np.clip(-np.tan(phi) * np.tan(delta), -1.0, 1.0))
    ra_mj = (24 * 60 / np.pi) * 0.0820 * dr * (
        ws * np.sin(phi) * np.sin(delta) + np.cos(phi) * np.cos(delta) * np.sin(ws)
    )
    ra_mm = 0.408 * ra_mj
    dt = np.clip(np.asarray(tmax, float) - np.asarray(tmin, float), 0.0, None)
    return np.clip(0.0023 * ra_mm * (np.asarray(tmean, float) + 17.8) * np.sqrt(dt), 0.0, None)


# --------------------------------------------------------------------------- #
# Plotting
# --------------------------------------------------------------------------- #
MONTH_STARTS = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
MONTH_LABELS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]


def panel(ax, cum_base, cum_now, n_days, title, ylabel, unit="mm"):
    x = np.arange(1, 366)
    p10, p25, p50, p75, p90 = [cum_base.quantile(q, axis=1) for q in (.1, .25, .5, .75, .9)]

    ax.fill_between(x, p10, p90, color="#8b949e", alpha=0.20, linewidth=0,
                    label=f"{BASELINE[0]}-{BASELINE[1]} 10th-90th pct.")
    ax.fill_between(x, p25, p75, color="#8b949e", alpha=0.32, linewidth=0,
                    label=f"{BASELINE[0]}-{BASELINE[1]} 25th-75th pct.")
    ax.plot(x, p50, color="white", ls="--", lw=2.0,
            label=f"{BASELINE[0]}-{BASELINE[1]} median")
    ax.plot(x[:n_days], cum_now[:n_days], color="#ff3b30", lw=2.6,
            label=f"{CURRENT_YEAR} (to DOY {n_days})")
    ax.plot([n_days], [cum_now[n_days - 1]], "o", color="#ff3b30", ms=5)

    ax.set_title("\n".join(textwrap.wrap(title, width=40)),
                 fontsize=12, fontweight="bold", color="white", pad=10)
    ax.set_ylabel(f"{ylabel} [{unit}]", fontsize=10, color="white")
    ax.set_xlim(1, 365)
    ax.set_xticks(MONTH_STARTS)
    ax.set_xticklabels(MONTH_LABELS, color="white")
    ax.tick_params(colors="white", labelsize=9, labelbottom=True)
    ax.grid(True, ls=":", alpha=0.25, color="gray")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#484f58")
    ax.legend(frameon=True, facecolor="#161b22", edgecolor="none",
              fontsize=8.5, labelcolor="white", loc="upper left")
    ax.set_box_aspect(1)   # square plotting area, independent of label sizes


def annotate(ax, cum_base, cum_now, n_days, unit="mm", signed=False):
    """Annotate with the current value, the baseline median, and the rank.

    `signed=True` for quantities that can cross zero (P - ET0). A "% of median"
    is only meaningful for a strictly positive quantity: -84 against a median of
    +324 is NOT "-26% of normal", it is 408 mm below normal. Reporting a ratio
    there invites a serious misreading, so signed quantities get an absolute
    difference instead.
    """
    now = float(cum_now[n_days - 1])
    med = float(cum_base.quantile(0.5, axis=1).iloc[n_days - 1])
    row = cum_base.iloc[n_days - 1].dropna()
    rank = int((row < now).sum()) + 1

    if signed or med <= 0 or now < 0:
        pct = float("nan")
        comp = f"{now - med:+.0f} {unit} vs. median"
    else:
        pct = 100 * now / med
        comp = f"{pct:.0f}% of median"

    txt = (f"{now:.0f} {unit}   |   median {med:.0f} {unit}   "
           f"({comp})   |   rank {rank}/{len(row) + 1}")
    ax.text(0.985, 0.04, txt, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=9, color="#c9d1d9",
            bbox=dict(facecolor="#161b22", edgecolor="#30363d", boxstyle="round,pad=0.4"))
    return now, med, pct, rank


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main(daily: pd.DataFrame | None = None, abbr: str | None = None,
         lat: float | None = None, station_label: str | None = None) -> dict:
    """Cumulative precipitation (and P - ET0) figure for the current year.

    With no arguments the series is downloaded from MeteoSwiss Open Data. Pass
    `daily`, `abbr` and `lat` to use a series already in memory — a cached copy,
    or the synthetic stand-in — through exactly the same computation.

    Returns:
        dict with end_doy, precip and water-balance totals, medians and ranks.
    """
    if daily is None:
        abbr, lat, _alt = resolve_station(STATION_NAME, STATION_ABBR)
        print("Downloading daily series (MeteoSwiss Open Data)...")
        daily = load_daily(abbr)
    label = station_label or STATION_NAME
    print(f"  {len(daily):,} daily records, {daily.index.min().date()} -> {daily.index.max().date()}")

    base_years = range(BASELINE[0], BASELINE[1] + 1)

    # --- precipitation -----------------------------------------------------
    p_base = doy_matrix(daily[P_PRECIP], base_years)
    p_now = doy_matrix(daily[P_PRECIP], range(CURRENT_YEAR, CURRENT_YEAR + 1))
    n_days = int(p_now[CURRENT_YEAR].last_valid_index() or 0)
    if n_days == 0:
        raise SystemExit(f"No {CURRENT_YEAR} data at station {abbr!r}.")
    print(f"  current year complete through aligned DOY {n_days}")
    report_gaps(p_base, "baseline precipitation")
    report_gaps(p_now, f"{CURRENT_YEAR} precipitation", upto=n_days)

    cum_p_base = cumulative(p_base)
    cum_p_now = cumulative(p_now)[CURRENT_YEAR].to_numpy()

    # --- optional water balance -------------------------------------------
    have_wb = INCLUDE_WATER_BALANCE and all(
        c in daily.columns for c in (P_TMEAN, P_TMAX, P_TMIN)
    )
    if INCLUDE_WATER_BALANCE and not have_wb:
        print("  [!] temperature parameters unavailable - skipping water-balance panel")

    if have_wb:
        pet = pd.Series(
            hargreaves_pet(daily[P_TMEAN], daily[P_TMAX], daily[P_TMIN],
                           daily.index.dayofyear, lat),
            index=daily.index,
        )
        wb = daily[P_PRECIP] - pet
        wb_base = doy_matrix(wb, base_years)
        wb_now = doy_matrix(wb, range(CURRENT_YEAR, CURRENT_YEAR + 1))
        cum_wb_base = cumulative(wb_base)
        cum_wb_now = cumulative(wb_now)[CURRENT_YEAR].to_numpy()

    # --- plot ---------------------------------------------------------------
    plt.style.use("dark_background")
    nrows = 2 if have_wb else 1
    fig, axes = plt.subplots(nrows, 1, figsize=(6.8, 7.6 * nrows), dpi=200,
                             sharex=True, constrained_layout=True)
    axes = np.atleast_1d(axes)

    panel(axes[0], cum_p_base, cum_p_now, n_days,
          f"{label} - cumulative precipitation", "Cumulative precipitation")
    stats = annotate(axes[0], cum_p_base, cum_p_now, n_days)
    print(f"\nPrecipitation to DOY {n_days}: {stats[0]:.0f} mm "
          f"({stats[2]:.0f}% of the {BASELINE[0]}-{BASELINE[1]} median), "
          f"rank {stats[3]} of {BASELINE[1]-BASELINE[0]+2} IN THE COMPARISON WINDOW "
          f"(not the full station record)")

    if have_wb:
        panel(axes[1], cum_wb_base, cum_wb_now, n_days,
              f"{label} - cumulative water balance P-ET0 (Hargreaves)",
              "Cumulative P - ET0")
        s2 = annotate(axes[1], cum_wb_base, cum_wb_now, n_days, signed=True)
        print(f"Water balance to DOY {n_days}: {s2[0]:.0f} mm "
              f"(median {s2[1]:.0f} mm), rank {s2[3]}")

    # The credit is the x label (see trajectory.make_figure).
    axes[-1].set_xlabel("Data: MeteoSchweiz Open Data (SwissMetNet), CC-BY",
                        loc="left", fontsize=7.5, color="#6e7681")
    out = figures_dir("station") / f"precip_{abbr}_{CURRENT_YEAR}.png"
    fig.savefig(out, dpi=200, facecolor=fig.get_facecolor())
    print(f"  -> {out}")
    plt.show()
    res = {"end_doy": n_days, "precip_mm": stats[0], "precip_median_mm": stats[1],
           "precip_pct_of_median": stats[2], "precip_rank": stats[3],
           "n_years": BASELINE[1] - BASELINE[0] + 2, "file": str(out)}
    if have_wb:
        res.update(wb_mm=s2[0], wb_median_mm=s2[1], wb_rank=s2[3])
    return res



