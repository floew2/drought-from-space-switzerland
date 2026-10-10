"""Synthetic cadastre: nine annual deliveries of grassland parcels.

The real Nutzungsflächen cannot be published, so the mask pipeline runs here on
parcels drawn to match the real master set (`targets.GEOMETRY`) and arranged so
that the things that broke the real pipeline are present too:

* outlines are REDRAWN every year (vertex jitter), fields are SPLIT into two and
  neighbours MERGED into one between deliveries, and object IDs are not stable
  across years — which is why parcels are matched by area overlap, not identity;
* some parcels were not grassland in some year, some changed management
  category, and Sömmerung (code 930) exists only in the latest delivery;
* the three sites differ: valley grassland at Linthebene and in the Rheintal,
  summering pasture at Alpstein.

Parcels are rectangles. For a w x h rectangle with aspect r = w / h the shape
index (perimeter / sqrt(area)) is s = 2(sqrt(r) + 1/sqrt(r)). To hit a target s,
solve u + 1/u = s/2 for u >= 1 and set r = u**2. (Setting r = s**2 / 16, as an
earlier attempt did, gives shape indices near 4.2 instead of 5.7.)
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import optimize, stats

from . import targets as T

MIN_HA = 0.25                 # the master set's own lower bound
YEARS = list(range(2017, 2026))
CRS = "EPSG:2056"

# Which LNF codes stand in for each component. A component is a set of codes
# that share a geometry distribution; the four analysis categories are unions
# of components and are derived from config.yaml exactly as for real data.
COMPONENT_CODES = {
    "I": [613, 650],                      # intensive (mown)
    "Em": [611, 612, 634, 635, 694],      # extensive, mown
    "Eg": [693],                          # extensive, grazed
    "Mo": [621, 622, 623, 935],           # mown, neither intensive nor extensive
    "Go": [616, 617, 618, 660, 625, 933], # grazed, other
    "S": [930],                           # Sömmerung — latest delivery only
}
NON_GRASS_CODES = [601, 602]              # Kunstwiese: temporary ley, excluded

# Master-set parcel counts per component. Only the category totals are known
# (targets.CATEGORY_TOTALS); the split of 693 out of extensive and of "mown
# other" out of mown is not, and is fixed here so that every total matches:
#   extensive = Em + Eg, mown = I + Em + Mo, grazed = Eg + Go + S.
_E_G = 400
MASTER_COUNTS = {
    "I": 23_758,
    "Em": 4_348 - _E_G,
    "Eg": _E_G,
    "Mo": 28_112 - 23_758 - (4_348 - _E_G),
    "Go": 6_474 - _E_G - T.FUNNEL["exempt_930"],
    "S": T.FUNNEL["exempt_930"],
}

# Share of each component that sits in each site. Valley grassland is mostly
# intensive; summering pasture is all at Alpstein.
SITE_SHARE = {
    "linthebene": {"I": 0.40, "Em": 0.45, "Eg": 0.30, "Mo": 0.30, "Go": 0.20, "S": 0.0},
    "rheintal":   {"I": 0.55, "Em": 0.45, "Eg": 0.30, "Mo": 0.50, "Go": 0.30, "S": 0.0},
    "alpstein":   {"I": 0.05, "Em": 0.10, "Eg": 0.40, "Mo": 0.20, "Go": 0.50, "S": 1.0},
}
# Cluster centre (LV95) and rotation of each site's synthetic landscape. The
# centres sit inside the real site boxes so per-site windows work unchanged;
# the geometry inside is entirely invented.
SITE_FRAME = {
    "linthebene": ((2_722_000.0, 1_227_500.0), 21.0),
    "rheintal":   ((2_763_500.0, 1_247_000.0), -14.0),
    "alpstein":   ((2_744_000.0, 1_229_000.0), 33.0),
}
SITE_ELEV = {"linthebene": (410, 440), "rheintal": (400, 520), "alpstein": (900, 2000)}

# Real funnel, as shares of the 2025 grassland population.
_N25 = T.FUNNEL["grassland_2025"]
P_NOT_GRASS = T.FUNNEL["dropped_not_grassland"] / (_N25 - T.FUNNEL["exempt_930"])
P_CHANGED = T.FUNNEL["dropped_changed_category"] / (_N25 - T.FUNNEL["exempt_930"])

# Between-delivery edits, per parcel and year.
P_SPLIT = 0.05
P_MERGE = 0.03
JITTER_M = 0.5


# --------------------------------------------------------------------------- #
# Distributions
# --------------------------------------------------------------------------- #
def _z(x, mu, sigma):
    return (np.log(x) - mu) / sigma


def _tln_cdf(x, mu, sigma, b=np.inf, a=MIN_HA):
    """CDF of a lognormal truncated to [a, b]."""
    fa, fb = stats.norm.cdf(_z(a, mu, sigma)), stats.norm.cdf(_z(b, mu, sigma))
    return (stats.norm.cdf(_z(x, mu, sigma)) - fa) / (fb - fa)


def _tln_mean(mu, sigma, b=np.inf, a=MIN_HA):
    za, zb = _z(a, mu, sigma), _z(b, mu, sigma)
    num = stats.norm.cdf(zb - sigma) - stats.norm.cdf(za - sigma)
    return np.exp(mu + sigma**2 / 2) * num / (stats.norm.cdf(zb) - stats.norm.cdf(za))


def _tln_ppf(q, mu, sigma, b=np.inf, a=MIN_HA):
    fa, fb = stats.norm.cdf(_z(a, mu, sigma)), stats.norm.cdf(_z(b, mu, sigma))
    return np.exp(mu + sigma * stats.norm.ppf(fa + q * (fb - fa)))


def fit_area(median_ha: float, lt_05: float, mean_ha: float | None = None
             ) -> tuple[float, float, float]:
    """Truncated lognormal with this median, P(< 0.5 ha) and, optionally, mean.

    The lower bound is fixed at 0.25 ha. Median and the small-parcel share fix
    the lognormal; if the mean is given too, an upper bound is fitted as the
    third parameter — an untruncated tail overshoots the real mean by a third
    for intensive meadows (2.4 ha against 1.77 ha).

    Args:
        median_ha: Target median parcel area in hectares.
        lt_05: Target share of parcels below 0.5 ha.
        mean_ha: Target mean area in hectares, or None for no upper bound.

    Returns:
        (mu, sigma, upper bound in ha), the bound being inf when not fitted.

    Raises:
        RuntimeError: If no distribution meets the targets.
    """
    def res2(p):
        mu, s = p[0], np.exp(p[1])
        return [_tln_cdf(median_ha, mu, s) - 0.5, _tln_cdf(0.5, mu, s) - lt_05]

    def res3(p):
        mu, s, b = p[0], np.exp(p[1]), median_ha + np.exp(p[2])
        return [_tln_cdf(median_ha, mu, s, b) - 0.5, _tln_cdf(0.5, mu, s, b) - lt_05,
                np.log(_tln_mean(mu, s, b)) - np.log(mean_ha)]

    sol = optimize.fsolve(res2, [np.log(median_ha), 0.0], full_output=True)
    if sol[2] != 1 or max(abs(r) for r in res2(sol[0])) > 1e-6:
        raise RuntimeError(f"area fit failed for median {median_ha}, <0.5 {lt_05}")
    mu, ls = sol[0]
    if mean_ha is None or _tln_mean(mu, np.exp(ls)) <= mean_ha:
        return float(mu), float(np.exp(ls)), float("inf")
    sol = optimize.fsolve(res3, [mu, ls, np.log(20 * median_ha)], full_output=True)
    if sol[2] != 1 or max(abs(r) for r in res3(sol[0])) > 1e-6:
        raise RuntimeError(f"bounded area fit failed for median {median_ha}")
    mu, ls, lb = sol[0]
    return float(mu), float(np.exp(ls)), float(median_ha + np.exp(lb))


def fit_grazed(eg: tuple[float, float, float], n: dict, sigma_s: float = 1.6
               ) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """Fit the 'other grazed' and Sömmerung components so the grazed MIXTURE matches.

    Grazed is a mixture of small extensive pastures (Eg, fixed), ordinary pasture
    (Go) and large summering units (S). The targets are the mixture's median,
    share below 0.5 ha, and mean area (total ha / parcels).

    Args:
        eg: (mu, sigma, upper bound) of the Eg component.
        n: Master counts per component.
        sigma_s: Log-scale spread of the summering units, fixed.

    Returns:
        ((mu, sigma, inf) for Go, (mu, sigma, inf) for S).
    """
    tgt = T.GEOMETRY["grazed"]
    cnt, ha = T.CATEGORY_TOTALS["grazed"]
    w = np.array([n["Eg"], n["Go"], n["S"]], float)
    w /= w.sum()

    def mix_cdf(x, p):
        return (w[0] * _tln_cdf(x, *eg) + w[1] * _tln_cdf(x, p[0], p[1])
                + w[2] * _tln_cdf(x, p[2], sigma_s))

    def res(p):
        p = [p[0], np.exp(p[1]), p[2]]
        mean = (w[0] * _tln_mean(*eg) + w[1] * _tln_mean(p[0], p[1])
                + w[2] * _tln_mean(p[2], sigma_s))
        return [mix_cdf(tgt["median_ha"], p) - 0.5,
                mix_cdf(0.5, p) - tgt["lt_05"],
                np.log(mean) - np.log(ha / cnt)]

    sol = optimize.fsolve(res, [np.log(0.8), np.log(0.9), np.log(4.0)], full_output=True)
    if sol[2] != 1 or max(abs(r) for r in res(sol[0])) > 1e-6:
        raise RuntimeError("grazed mixture fit failed")
    mu_g, ls_g, mu_s = sol[0]
    inf = float("inf")
    return (float(mu_g), float(np.exp(ls_g)), inf), (float(mu_s), sigma_s, inf)


def shape_sigma(median_shape: float) -> float:
    """Log-spread of (shape - 4) such that P(shape > 16) matches the real tail.

    The shape index is modelled as 4 + X, X lognormal with median
    (median_shape - 4). Its spread is set by the one tail statistic known from
    the real data: about two parcels in 5,000 exceed 16.
    """
    z = stats.norm.isf(T.SHAPE_TAIL_GT16)
    return float(np.log(12.0 / (median_shape - 4.0)) / z)


def aspect_for_shape(s: np.ndarray) -> np.ndarray:
    """Aspect ratio r = w/h (>= 1) of a rectangle whose shape index is s.

    s = 2(sqrt(r) + 1/sqrt(r)). With u = sqrt(r): u + 1/u = s/2, u >= 1.
    """
    half = np.asarray(s, float) / 2.0
    u = (half + np.sqrt(np.maximum(half**2 - 4.0, 0.0))) / 2.0
    return u**2


@dataclass(frozen=True)
class Component:
    """Geometry distribution of one component: area and shape index."""
    area: tuple[float, float, float]   # (mu, sigma, upper bound) of the truncated lognormal, ha
    shape_median: float
    shape_sigma: float


def components() -> dict[str, Component]:
    """Calibrated geometry distribution for every component."""
    g = T.GEOMETRY
    ci, ce = T.CATEGORY_TOTALS["intensive"], T.CATEGORY_TOTALS["extensive"]
    ai = fit_area(g["intensive"]["median_ha"], g["intensive"]["lt_05"], ci[1] / ci[0])
    ae = fit_area(g["extensive"]["median_ha"], g["extensive"]["lt_05"], ce[1] / ce[0])
    go, s = fit_grazed(ae, MASTER_COUNTS)
    si, se, sg = (g[k]["median_shape"] for k in ("intensive", "extensive", "grazed"))
    return {
        "I": Component(ai, si, shape_sigma(si)),
        "Em": Component(ae, se, shape_sigma(se)),
        "Eg": Component(ae, se, shape_sigma(se)),
        "Mo": Component(ai, si, shape_sigma(si)),
        "Go": Component(go, sg, shape_sigma(sg)),
        "S": Component(s, sg, shape_sigma(sg)),
    }


# --------------------------------------------------------------------------- #
# The 2025 population and its layout
# --------------------------------------------------------------------------- #
def draw_population(rng: np.random.Generator, scale: float) -> pd.DataFrame:
    """Draw the 2025 grassland parcels: component, code, site, status, w, h.

    Counts are the master counts scaled, inflated by the real drop rate so that
    after the persistence test the MASTER set has the target counts.
    """
    comp = components()
    keep = 1.0 - P_NOT_GRASS - P_CHANGED
    rows = []
    for c, n_master in MASTER_COUNTS.items():
        n = round(n_master * scale / (1.0 if c == "S" else keep))
        # Stratified quantiles, then shuffled: the sample hits the target
        # distribution closely even at a few hundred parcels.
        q = (rng.permutation(n) + rng.uniform(0, 1, n)) / n
        area_ha = _tln_ppf(q, *comp[c].area)
        x = np.exp(np.log(comp[c].shape_median - 4.0)
                   + comp[c].shape_sigma * stats.norm.ppf(
                       (rng.permutation(n) + rng.uniform(0, 1, n)) / n))
        shape = 4.0 + x
        r = aspect_for_shape(shape)
        h = np.sqrt(area_ha * 1e4 / r)
        sites = rng.choice(list(SITE_SHARE), size=n,
                           p=_norm([SITE_SHARE[s][c] for s in SITE_SHARE]))
        for i in range(n):
            rows.append({"component": c, "code": int(rng.choice(COMPONENT_CODES[c])),
                         "site": sites[i], "w": h[i] * r[i], "h": h[i]})
    df = pd.DataFrame(rows)

    status = np.full(len(df), "stable", dtype=object)
    s930 = df.component.eq("S").to_numpy()
    status[s930] = "exempt_930"
    u = rng.uniform(size=len(df))
    status[~s930 & (u < P_NOT_GRASS)] = "not_grassland"
    status[~s930 & (u >= P_NOT_GRASS) & (u < P_NOT_GRASS + P_CHANGED)] = "changed_category"
    df["status"] = status
    return df


def _norm(p):
    p = np.asarray(p, float)
    return p / p.sum()


def layout(df: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    """Place every parcel in its site's landscape: rows of fields, no overlaps.

    Shelf packing: parcels are binned by depth (short side) and shuffled within
    a bin, so similar-sized fields share a row while management categories mix
    — an extensive meadow usually borders intensive ones, as it does on the
    ground. Fields in a row abut with probability 0.6 (a shared boundary) and are
    otherwise separated by a track or hedge of 3-12 m.

    Adds columns x, y (lower-left corner, local metres), right_nb (index of an
    abutting right-hand neighbour or -1), and elev_m.
    """
    df = df.copy()
    df["x"] = np.nan
    df["y"] = np.nan
    df["right_nb"] = -1
    for site in SITE_SHARE:
        idx = df.index[df.site == site].to_numpy()
        if not len(idx):
            continue
        sub = df.loc[idx]
        width = max(np.sqrt(sub.w.mul(sub.h).sum() / 0.55), sub.w.max() * 1.05)
        bins = pd.qcut(sub.h.rank(method="first"), min(25, len(sub)), labels=False)
        order = []
        for b in sorted(bins.unique(), reverse=True):
            members = sub.index[bins == b].to_numpy()
            order.extend(rng.permutation(members))
        x = y = 0.0
        row_h = 0.0
        prev, abuts = -1, False
        for i in order:
            w, h = df.at[i, "w"], df.at[i, "h"]
            if x + w > width and x > 0:
                y += row_h + (0.0 if rng.uniform() < 0.5 else rng.uniform(3, 10))
                x, row_h, prev = 0.0, 0.0, -1
            df.at[i, "x"], df.at[i, "y"] = x, y
            row_h = max(row_h, h)
            if prev >= 0 and abuts:
                df.at[prev, "right_nb"] = i
            abuts = rng.uniform() < 0.6
            x += w + (0.0 if abuts else rng.uniform(3, 12))
            prev = i
        lo, hi = SITE_ELEV[site]
        e = rng.uniform(lo, hi, len(idx))
        df.loc[idx, "elev_m"] = np.where(df.loc[idx, "component"].eq("S"),
                                         np.maximum(e, 1400), e).round(0)
        # centre the cluster on its frame origin
        cx = df.loc[idx, "x"].add(df.loc[idx, "w"]).max() / 2
        cy = df.loc[idx, "y"].add(df.loc[idx, "h"]).max() / 2
        df.loc[idx, "x"] -= cx
        df.loc[idx, "y"] -= cy
    return df


# --------------------------------------------------------------------------- #
# Annual deliveries
# --------------------------------------------------------------------------- #
def _codes_by_year(df: pd.DataFrame, rng: np.random.Generator) -> np.ndarray:
    """Code per parcel and year (len(YEARS) columns); 0 = absent from the delivery."""
    codes = np.repeat(df["code"].to_numpy()[:, None], len(YEARS), axis=1)
    past = len(YEARS) - 1
    other_cat = {"I": "Em", "Em": "I", "Eg": "Go", "Mo": "Go", "Go": "Mo"}
    for k, (st, comp) in enumerate(zip(df.status, df.component)):
        if st == "exempt_930":
            codes[k, :past] = 0
        elif st in ("not_grassland", "changed_category"):
            yrs = rng.choice(past, size=int(rng.integers(1, 4)), replace=False)
            for y in yrs:
                codes[k, y] = (int(rng.choice(NON_GRASS_CODES)) if st == "not_grassland"
                               else int(rng.choice(COMPONENT_CODES[other_cat[comp]])))
    return codes


def deliveries(df: pd.DataFrame, rng: np.random.Generator):
    """Build the nine annual deliveries as GeoDataFrames in EPSG:2056.

    Args:
        df: Laid-out population (from `layout`).
        rng: Seeded generator.

    Returns:
        dict year -> GeoDataFrame with columns objectid, lnf_code, geometry; and
        the 2025 frame's objectid per parcel (the identity the cut series use).
    """
    import geopandas as gpd
    from shapely import affinity
    from shapely.geometry import Polygon, box

    codes = _codes_by_year(df, rng)
    out, ids_2025 = {}, None
    for yi, year in enumerate(YEARS):
        latest = year == YEARS[-1]
        absorbed = {}                              # right neighbour -> owner
        if not latest:
            used = set()
            for k, nb in enumerate(df.right_nb.to_numpy()):
                if (nb >= 0 and codes[k, yi] and codes[k, yi] == codes[nb, yi]
                        and k not in used and nb not in used
                        and rng.uniform() < P_MERGE):
                    absorbed[nb] = k
                    used.update((k, nb))
        partner_of = {o: m for m, o in absorbed.items()}
        items = []                                 # (site, polygon, parcel)
        for k, r in enumerate(df.itertuples()):
            if not codes[k, yi] or k in absorbed:
                continue
            if k in partner_of:
                # Merged with its abutting right-hand neighbour into one field.
                # Built from vertices, not by union: the shared edge is exact
                # only up to float rounding, and a union would leave a sliver.
                p = df.iloc[partner_of[k]]
                x1 = p.x + p.w
                items.append((r.site, Polygon([(r.x, r.y), (x1, r.y), (x1, r.y + p.h),
                                      (p.x, r.y + p.h), (p.x, r.y + r.h),
                                      (r.x, r.y + r.h)]).buffer(0), k))
            elif not latest and rng.uniform() < P_SPLIT:
                # Any field can be split, small ones included; a part is never
                # below 30 % of the field, so never below the 0.05 ha sliver
                # floor. (Under the old 0.25 ha filter on past deliveries, such
                # splits wrongly failed their parent — METHOD §6, #11.)
                t = rng.uniform(0.3, 0.7) * r.w
                items.append((r.site, box(r.x, r.y, r.x + t, r.y + r.h), k))
                items.append((r.site, box(r.x + t, r.y, r.x + r.w, r.y + r.h), k))
            else:
                items.append((r.site, box(r.x, r.y, r.x + r.w, r.y + r.h), k))
        polys, gcodes, owner = [], [], []
        for site, g, k in items:
            (cx, cy), ang = SITE_FRAME[site]
            g = affinity.rotate(_jitter(g, rng), ang, origin=(0, 0))
            polys.append(affinity.translate(g, cx, cy))
            gcodes.append(int(codes[k, yi]))
            owner.append(k)
        # Object IDs are reassigned every year, as in the real deliveries:
        # nothing that follows an ID through time can work.
        oid = rng.permutation(len(polys)) + 100_000 + 10_000 * yi
        gdf = gpd.GeoDataFrame({"objectid": oid, "lnf_code": gcodes},
                               geometry=polys, crs=CRS)
        out[year] = gdf
        if latest:
            ids_2025 = pd.Series(oid, index=owner)
    return out, ids_2025


def _jitter(geom, rng):
    """Move every vertex a little: the same field, redrawn."""
    from shapely.geometry import Polygon

    ring = np.asarray(geom.exterior.coords)[:-1]
    ring = ring + rng.normal(0, JITTER_M, ring.shape)
    out = Polygon(ring)
    return out if out.is_valid else geom


def generate(rng: np.random.Generator, scale: float):
    """Population, layout and deliveries in one call.

    Returns:
        (parcels, deliveries): `parcels` has one row per 2025 parcel with its
        truth (component, site, status, elevation, objectid), `deliveries`
        maps year -> GeoDataFrame.
    """
    pop = layout(draw_population(rng, scale), rng)
    dl, ids = deliveries(pop, rng)
    pop["objectid"] = ids.reindex(pop.index).to_numpy()
    pop["area_ha"] = (pop.w * pop.h / 1e4).round(4)
    return pop, dl
