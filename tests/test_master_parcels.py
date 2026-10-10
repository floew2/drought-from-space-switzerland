"""The master parcel set: overlap matching and the sliver floor (METHOD §6, #11)."""

from __future__ import annotations

import geopandas as gpd
from shapely.geometry import box

from drought import mask


def _deliveries(tmp_path, years):
    for year, geoms in years.items():
        d = tmp_path / str(year)
        d.mkdir()
        gpd.GeoDataFrame({"lnf_code": [613] * len(geoms)}, geometry=geoms,
                         crs="EPSG:2056").to_file(d / f"nutzungen_{year}.shp")
    return tmp_path


def test_redrawn_or_split_small_field_keeps_its_history(tmp_path):
    x0, y0 = 2_722_000, 1_227_500
    field_2025 = box(x0, y0, x0 + 130, y0 + 20)              # 0.26 ha
    redrawn = box(x0, y0, x0 + 120, y0 + 20)                 # 0.24 ha, same ground
    split = [box(x0, y0, x0 + 60, y0 + 20), box(x0 + 60, y0, x0 + 130, y0 + 20)]
    shp = _deliveries(tmp_path, {2023: split, 2024: [redrawn], 2025: [field_2025]})

    kept = mask.master_parcels(shp, verbose=False)
    assert len(kept) == 1, "a stable field must survive a redraw or split below 0.25 ha"

    # The old behaviour, kept reproducible: a 0.25 ha floor on every delivery
    # loses the same field.
    old = mask.master_parcels(shp, verbose=False, sliver_floor_m2=2500)
    assert len(old) == 0


def test_slivers_are_dropped(tmp_path):
    x0, y0 = 2_722_000, 1_227_500
    sliver = box(x0, y0, x0 + 20, y0 + 20)                   # 0.04 ha
    shp = _deliveries(tmp_path, {2024: [sliver], 2025: [sliver]})
    assert len(mask.master_parcels(shp, verbose=False)) == 0
