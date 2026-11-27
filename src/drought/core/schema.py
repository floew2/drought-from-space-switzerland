"""The shape of config.yaml, checked when it is loaded.

Every section is a model that forbids unknown keys, so a misspelt or misplaced
setting fails at load time, naming the key, instead of being silently ignored.
The models only check: the configuration is still used as the plain attribute
tree in `config.CFG`, so overriding a value for one run stays simple.
"""

from __future__ import annotations

from pathlib import Path
from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field, ValidationError, field_validator

__all__ = ["Settings", "validate"]

Doy = Annotated[int, Field(ge=1, le=366)]
DoyWindow = tuple[Doy, Doy]
Share = Annotated[float, Field(ge=0, le=1)]
Percent = Annotated[float, Field(ge=0, le=100)]
Index = Literal["NDVI", "kNDVI", "NDMI", "NDRE"]


class _Section(BaseModel):
    model_config = ConfigDict(extra="forbid")


class Project(_Section):
    current_year: int
    baseline_years: list[int] = Field(min_length=3)


class Paths(_Section):
    scenes: str
    series: str
    labels: str
    masks: str
    parcels: str
    parcel_ids: str
    station: str
    figures: str
    tables: str


class Site(_Section):
    label: str
    note: str = ""
    bbox_wgs84: tuple[float, float, float, float]
    stratum: Literal["valley", "summering"] = "valley"

    @field_validator("bbox_wgs84")
    @classmethod
    def _inside_switzerland(cls, v):
        w, s, e, n = v
        if not (5.9 <= w < e <= 10.6 and 45.8 <= s < n <= 47.9):
            raise ValueError("bbox must be (west, south, east, north) inside Switzerland: "
                             "swissEO covers Switzerland only")
        return v


class Landuse(_Section):
    mown: list[int]
    grazed: list[int]
    intensive: list[int]
    extensive: list[int]
    earliest_legal_cut_doy: Doy


class Mask(_Section):
    res: float = Field(gt=0)
    sliver_floor_m2: float = Field(ge=0)
    subpixel: int = Field(ge=1)
    min_overlap: Share
    min_years_same: int = Field(ge=0)
    simplify_m: float = Field(ge=0)
    psf_mtf_nyquist: float = Field(gt=0, lt=1)
    psf_model: Literal["box_gauss", "gauss"]
    clip_bbox_lv95: tuple[float, float, float, float] | None
    purity_min: Share


class Clean(_Section):
    min_parcel_coverage: Share
    haze_window: int = Field(ge=2)
    haze_ndvi_z: float
    haze_gap: float


class Detect(_Section):
    index: Index
    confirm: Index
    season: DoyWindow
    truncate_to_current: bool
    min_drop_abs: float
    min_drop_rate: float
    drop_fraction: float
    lookback_days: int
    rise_window: int
    min_rise_rate: float
    merge_days: int
    persist_frac: float
    persist_always: bool
    min_concentration: float
    require_followup: int
    confirm_frac: float
    amp_window: int


class Changepoint(_Section):
    penalty_k: float
    min_seg_days: int
    min_jump: float


class Analysis(_Section):
    merge_regimes: bool
    regimes: list[Literal["intensive", "extensive", "mown", "grazed"]]
    late_split_doy: Doy
    integration_window: DoyWindow
    min_window_coverage: Share
    min_obs_per_year: int = Field(ge=2)
    max_gap_days: int = Field(ge=1)
    smooth_days: int = Field(ge=0)
    min_coverage: Share


class Evaluate(_Section):
    match_tolerance_days: int = Field(ge=0)


class Precip(_Section):
    station_name: str
    baseline: tuple[int, int]
    include_water_balance: bool


class Usable(_Section):
    cloud_clear: list[Literal[0, 2]]
    terrain_drop: list[int]
    terrain_max_incidence: float = Field(gt=0, le=180)
    scl_drop: list[int]


class Vhi(_Section):
    collection: str
    asset: str
    nodata: list[int]


class Rgb(_Section):
    year: int
    doy: DoyWindow


class Swisseo(_Section):
    stac: str
    collection: str
    scale: float
    offset: float
    season_doy: DoyWindow
    workers: int = Field(ge=1, le=16)
    usable: Usable
    vhi: Vhi
    rgb: Rgb


class Parcels(_Section):
    year: int


class Condition(_Section):
    doy: DoyWindow
    sites: list[str]
    min_pixels: int = Field(ge=1)


class Sensitivity(_Section):
    index: Index
    n_classes: int = Field(ge=2)
    strict_purity: Share
    max_site_cloud: Percent


class Cuts(_Section):
    sites: list[str]
    doy: DoyWindow
    sample_n: int = Field(ge=1)
    min_area_ha: float = Field(ge=0)
    edge_buffer_m: float = Field(le=0)
    seed: int
    max_site_cloud: Percent


class Anomaly(_Section):
    index: Index
    sites: list[str]
    min_obs: int = Field(ge=1)
    min_baseline_years: int = Field(ge=2)


class Synthetic(_Section):
    seed: int
    scale: float = Field(gt=0, le=1)
    paths: Paths


class Provenance(_Section):
    require: bool


class Settings(_Section):
    """The whole of config.yaml."""
    project: Project
    paths: Paths
    sites: dict[str, Site]
    windows: dict[str, DoyWindow]
    landuse: Landuse
    mask: Mask
    clean: Clean
    detect: Detect
    changepoint: Changepoint
    analysis: Analysis
    evaluate: Evaluate
    precip: Precip
    swisseo: Swisseo
    parcels: Parcels
    condition: Condition
    sensitivity: Sensitivity
    cuts: Cuts
    anomaly: Anomaly
    synthetic: Synthetic
    provenance: Provenance

    @field_validator("windows")
    @classmethod
    def _has_default(cls, v):
        if "default" not in v:
            raise ValueError("windows needs a `default` entry")
        return v


def validate(raw: dict, source: str | Path = "config.yaml") -> Settings:
    """Check a loaded configuration against the schema.

    Raises:
        ValueError: Listing every key that is missing, unknown or out of range.
    """
    try:
        s = Settings.model_validate(raw)
    except ValidationError as exc:
        lines = [f"  {'.'.join(map(str, e['loc']))}: {e['msg']}" for e in exc.errors()]
        raise ValueError(f"{Path(source).name} does not match the schema:\n"
                         + "\n".join(lines)) from None
    known = set(s.sites)
    for section in ("condition", "cuts", "anomaly"):
        bad = [x for x in getattr(s, section).sites if x not in known]
        if bad:
            raise ValueError(f"{Path(source).name}: {section}.sites names unknown sites {bad}")
    bad = [k for k in s.windows if k != "default" and k not in known]
    if bad:
        raise ValueError(f"{Path(source).name}: windows names unknown sites {bad}")
    return s
