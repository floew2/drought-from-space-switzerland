"""Figures for the grassland mask.

Kept apart from the mask computation so that building the mask never needs
matplotlib, and a figure never rebuilds the parcel set to compute numbers it
could be handed (METHOD §6, bug 5).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .config import mask_path
from .mask import GROUPS, PURITY_MIN, RES, master_parcels, psf_caption, psf_params
from .mask.build import area_table, verify_mask
from .mask.io import site_bounds_lv95
from .mask.purity import _burn

__all__ = ["area_check", "figure_mask_overview", "figure_parcel_size", "figure_psf",
           "figure_purity_explainer", "figure_purity_map"]


def figure_purity_map(site: str, group: str = "mown", out: str | Path | None = None,
                      inset: tuple[float, float, float, float] | None = None,
                      inset_frac: float = 0.14, threshold: float | None = None) -> np.ndarray | None:
    """Continuous purity for one landscape, with a zoomed inset.

    Deliberately not thresholded: the gradient is the point. The threshold is
    only a marker on the colour bar, so a reader sees where it sits in the
    distribution. Real-data maps are not publishable (METHOD §7, gallery rules).

    Args:
        site: Site key in config.yaml.
        group: Mask band: "mown", "grazed", "intensive" or "extensive".
        out: Output PNG (default: figures/mask/purity_map_<site>_<group>.png).
        inset: (x0, y0, x1, y1) in LV95; default is the densest window.
        inset_frac: Inset size as a share of the extent, when chosen automatically.
        threshold: Threshold to mark on the colour bar (default: config).

    Returns:
        The cropped purity array, or None if the group is absent at the site.
    """
    import matplotlib.pyplot as plt
    import rasterio
    from matplotlib.patches import Rectangle

    from .config import figures_dir

    threshold = PURITY_MIN if threshold is None else threshold
    band = list(GROUPS).index(group) + 1
    x0, y0, x1, y1 = site_bounds_lv95(site)
    verify_mask()
    with rasterio.open(mask_path()) as ds:
        win = rasterio.windows.from_bounds(x0, y0, x1, y1, ds.transform)
        a = ds.read(band, window=win)
        tr = ds.window_transform(win)
    if not a.size or not (a > 0).any():
        print(f"no {group} parcels in {site}")
        return None

    res = tr.a
    # Crop to where this group actually occurs, otherwise a landscape whose
    # parcels sit in one corner renders as mostly empty black.
    rows, cols = np.where(a > 0)
    pad = 8
    r0, r1 = max(rows.min() - pad, 0), min(rows.max() + pad + 1, a.shape[0])
    c0, c1 = max(cols.min() - pad, 0), min(cols.max() + pad + 1, a.shape[1])
    a = a[r0:r1, c0:c1]
    ox, oy = tr.c + c0 * res, tr.f - r0 * res
    ext = [ox, ox + a.shape[1] * res, oy - a.shape[0] * res, oy]

    # Densest window, so the inset shows fields rather than empty ground.
    if inset is None:
        n = max(int(min(a.shape) * inset_frac), 12)
        occ = (a > 0).astype(float)
        cs = occ.cumsum(0).cumsum(1)
        best, br, bc = -1, 0, 0
        for r in range(0, a.shape[0] - n, max(n // 3, 1)):
            for c in range(0, a.shape[1] - n, max(n // 3, 1)):
                tot = (cs[r + n, c + n] - cs[r, c + n] - cs[r + n, c] + cs[r, c]
                       if r and c else occ[r:r + n, c:c + n].sum())
                if tot > best:
                    best, br, bc = tot, r, c
        inset = (ext[0] + bc * res, ext[3] - (br + n) * res,
                 ext[0] + (bc + n) * res, ext[3] - br * res)
    ir0 = int((ext[3] - inset[3]) / res); ic0 = int((inset[0] - ext[0]) / res)
    ir1 = int((ext[3] - inset[1]) / res); ic1 = int((inset[2] - ext[0]) / res)
    sub = a[ir0:ir1, ic0:ic1]

    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(11, 8.4), dpi=180, constrained_layout=True)
    m = np.ma.masked_where(a <= 0, a)
    # Black -> white, so purity reads as brightness: pure cells glow, mixed
    # edges fade out. A hue ramp made the rim a colour rather than a gradient.
    cmap = plt.cm.gray.copy(); cmap.set_bad("black")
    im = ax.imshow(m, extent=ext, cmap=cmap, vmin=0, vmax=1,
                   interpolation="nearest")
    ax.add_patch(Rectangle((inset[0], inset[1]), inset[2] - inset[0],
                           inset[3] - inset[1], fill=False, edgecolor="#58a6ff",
                           lw=1.8, zorder=5))

    # --- inset, upper right
    axi = ax.inset_axes([0.615, 0.615, 0.375, 0.375])
    mi = np.ma.masked_where(sub <= 0, sub)
    axi.imshow(mi, cmap=cmap, vmin=0, vmax=1, interpolation="nearest")
    axi.set_xticks([]); axi.set_yticks([])
    for sp in axi.spines.values():
        sp.set_color("#58a6ff"); sp.set_linewidth(1.8)
    axi.text(0.02, 0.975, f"{(inset[2]-inset[0])/1000:.1f} km wide",
             transform=axi.transAxes, fontsize=8, color="#58a6ff", va="top",
             bbox=dict(facecolor="black", edgecolor="none", pad=2))

    # --- scale bar, lower left
    span = ext[1] - ext[0]
    bar = max(round(span / 5 / 1000), 1) * 1000
    bx, by = ext[0] + span * 0.04, ext[2] + (ext[3] - ext[2]) * 0.055
    ax.plot([bx, bx + bar], [by, by], color="white", lw=3, solid_capstyle="butt",
            zorder=6, path_effects=[])
    ax.text(bx + bar / 2, by + (ext[3] - ext[2]) * 0.018, f"{bar / 1000:g} km",
            ha="center", fontsize=9, color="white", zorder=6,
            bbox=dict(facecolor="black", edgecolor="none", pad=1.5))

    # --- north arrow, lower right
    nx, ny = ext[1] - span * 0.055, ext[2] + (ext[3] - ext[2]) * 0.055
    ax.annotate("", xy=(nx, ny + (ext[3] - ext[2]) * 0.075), xytext=(nx, ny),
                zorder=6,
                arrowprops=dict(arrowstyle="-|>", color="white", lw=1.6))
    ax.text(nx, ny - (ext[3] - ext[2]) * 0.012, "N", ha="center", va="top",
            fontsize=11, color="white", fontweight="bold", zorder=6,
            bbox=dict(facecolor="black", edgecolor="none", pad=1.5))

    # --- LV95 coordinates, sparse
    ax.set_xlabel("E [m], CH1903+ / LV95", fontsize=9, color="#8b949e")
    ax.set_ylabel("N [m]", fontsize=9, color="#8b949e")
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    ax.ticklabel_format(style="plain", useOffset=False)
    ax.tick_params(colors="#8b949e", labelsize=8)
    for sp in ax.spines.values():
        sp.set_color("#484f58")

    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.015)
    cb.set_label(f"purity — share of the {RES:g} m cell's response from inside a parcel",
                 fontsize=9, color="white")
    cb.ax.axhline(threshold, color="#ff3b30", lw=2)
    cb.ax.tick_params(colors="white", labelsize=8)
    cb.ax.text(1.6, threshold, f" {threshold:.2f}", color="#ff3b30",
               fontsize=8.5, va="center", transform=cb.ax.get_yaxis_transform())

    kept = 100 * (a[a > 0] >= threshold).mean()
    ax.set_title(f"Pixel purity  ·  {site}  ·  {group}", fontsize=12,
                 fontweight="bold", color="#ff3b30", loc="left")
    fig.supxlabel(f"Continuous, not thresholded — the gradient is the point. "
                  f"Black = no parcel of this group. At {threshold:.2f}, "
                  f"{kept:.0f}% of cells are retained.",
                  fontsize=7.5, color="#6e7681")
    out = out or (figures_dir("mask") / f"purity_map_{site}_{group}.png")
    fig.savefig(out, dpi=180, facecolor=fig.get_facecolor())
    print(f"  -> {out}")
    plt.show()
    return a


def area_check(site: str | None = None,
               thresholds: tuple[float, ...] = (0.0, 0.25, 0.5, 0.75, 0.9),
               out: str | Path | None = None, master=None) -> dict:
    """Area cost of the purity threshold: the table, and the curve per category.

    Args:
        site: Restrict to a site window, or None for the full raster.
        thresholds: Thresholds for the printed table.
        out: Output PNG (default: figures/mask/area_check_<site|all>.png).
        master: Master parcel set; built if None (slow).

    Returns:
        The `area_table` dict, without the plotting window.
    """
    import matplotlib.pyplot as plt
    import rasterio

    from .config import figures_dir

    out_d = area_table(master, site, thresholds)
    win = out_d.pop("_window")
    names = list(GROUPS)
    cell_ha = (RES * RES) / 1e4

    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(8.0, 5.6), dpi=180, constrained_layout=True)
    col = {"intensive": "#3fb950", "extensive": "#58a6ff",
           "mown": "#d29922", "grazed": "#8250df"}
    grid_t = np.round(np.arange(0.0, 0.96, 0.05), 2)
    verify_mask()
    with rasterio.open(mask_path()) as ds:
        for i, n in enumerate(names, 1):
            a = ds.read(i, window=win)
            poly = out_d[n]["polygon_ha"]
            if not poly:
                continue
            r = [float((a >= max(t, 1e-9)).sum()) * cell_ha / poly
                 for t in grid_t]
            ax.plot(grid_t, r, "-o", ms=3.5, lw=2.2, color=col[n],
                    label=f"{n}  ({out_d[n]['median_ha']:.2f} ha, "
                          f"shape {out_d[n]['median_shape']:.1f})")
    ax.axhline(1.0, color="#8b949e", ls=":", lw=1.4)
    ax.axvline(PURITY_MIN, color="#ff3b30", ls="--", lw=1.8,
               label=f"threshold in use ({PURITY_MIN:.2f})")
    ax.set_xlabel("purity threshold", fontsize=10, color="white")
    ax.set_ylabel("observed area / polygon area", fontsize=10, color="white")
    ax.set_xlim(0, 0.97)
    ax.set_title(f"How much of the master parcel set do we measure?"
                 f"  ·  {site or 'full extent'}", fontsize=11.5,
                 fontweight="bold", color="#ff3b30", loc="left")
    ax.legend(fontsize=8.5, frameon=True, facecolor="#161b22",
              edgecolor="none", labelcolor="white", loc="upper right")
    ax.grid(alpha=0.22, ls=":", color="gray")
    ax.tick_params(colors="white", labelsize=9)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color("#484f58")
    fig.supxlabel("Above 1 at low thresholds is the PSF halo, and is required. "
                  "Crossing 1 is not the target — it happens where about half "
                  "the retained\ncells are contaminated. Legend gives median "
                  "parcel size and shape index (4.0 = a perfect square).",
                  fontsize=7, color="#6e7681")
    out = out or (figures_dir("mask") / f"area_check_{site or 'all'}.png")
    fig.savefig(out, dpi=180, facecolor=fig.get_facecolor())
    print(f"\n  -> {out}")
    plt.show()
    return out_d


def figure_purity_explainer(site: str = "linthebene", group: str = "intensive",
                            out: str | Path | None = None, zoom_m: float = 400,
                            thresholds: list[float] | None = None, master=None,
                            areas: dict[str, float] | None = None,
                            zoom_levels: tuple[float, ...] = (0.0, 0.6, 1.0)) -> dict:
    """Why purity exists and what it costs, in two panels.

    A: the same patch at several thresholds, so the fields visibly erode inward.
    B: observed / mapped area against the threshold, per category. Pass `areas`
    (the four polygon totals `build_mask` prints) so the figure never rebuilds the
    parcel set for four numbers (METHOD §6, bug 5).

    Args:
        site: Site key for panel A.
        group: Mask band for panel A.
        out: Output PNG (default: figures/mask/purity_explainer_<site>_<group>.png).
        zoom_m: Width of the panel-A window in metres.
        thresholds: Thresholds for the panel-B curve.
        master: Master parcel set, used only if `areas` is not given.
        areas: Polygon area per category in ha.
        zoom_levels: Thresholds shown in panel A; the one equal to
            `purity_min` is framed in red.

    Returns:
        dict category -> list of area ratios over `thresholds`.
    """
    import matplotlib.pyplot as plt
    import rasterio
    from matplotlib.patches import Rectangle

    from .config import figures_dir

    thresholds = list(thresholds or np.round(np.arange(0.0, 0.96, 0.05), 2))
    names = list(GROUPS)

    # --- data -------------------------------------------------------------
    x0, y0, x1, y1 = site_bounds_lv95(site)
    verify_mask()
    with rasterio.open(mask_path()) as ds:
        b_ = ds.bounds
        xa, ya = max(b_.left, x0), max(b_.bottom, y0)
        xb, yb = min(b_.right, x1), min(b_.top, y1)
        band = names.index(group) + 1
        win = rasterio.windows.from_bounds(xa, ya, xb, yb, ds.transform)
        big = ds.read(band, window=win)
        full = {n: ds.read(names.index(n) + 1,
                           window=rasterio.windows.from_bounds(
                               b_.left, b_.bottom, b_.right, b_.top,
                               ds.transform))
                for n in names}

    # densest window for the zoom
    npx = max(int(zoom_m / RES), 20)
    occ = (big > 0).astype(float)
    best, br, bc = -1, 0, 0
    step = max(npx // 3, 1)
    for r_ in range(0, max(big.shape[0] - npx, 1), step):
        for c_ in range(0, max(big.shape[1] - npx, 1), step):
            v = occ[r_:r_ + npx, c_:c_ + npx].sum()
            if v > best:
                best, br, bc = v, r_, c_
    zoom = big[br:br + npx, bc:bc + npx]

    # Panel C needs ONE number per category: the mapped polygon area. Nothing
    # else here touches polygons — panel A is the Gaussian, panel B reads the
    # raster.
    #
    # Those four numbers are printed by main(). Pass them and the figure never
    # goes near a shapefile:
    #     areas={"mown": 44628, "grazed": 36223,
    #            "intensive": 42034, "extensive": 2585}
    # Or pass master= to sum them from a set you already built. Only with
    # neither does it rebuild, and rebuilding means eight yearly overlays —
    # minutes to hours for four numbers you already have.
    cell_ha = (RES * RES) / 1e4
    if areas is None:
        if master is None:
            print("  building the parcel set just to total four areas — pass "
                  "areas= or master= to skip this")
            master = master_parcels(verbose=False)
        areas = {n: float(master[master["_code"].isin(GROUPS[n])]
                          .geometry.area.sum() / 1e4) for n in names}

    curves = {}
    for n in names:
        poly = float(areas.get(n, 0.0))
        if poly <= 0:
            continue
        curves[n] = [float((full[n] >= max(t, 1e-9)).sum()) * cell_ha / poly
                     for t in thresholds]

    # --- figure -----------------------------------------------------------
    plt.style.use("dark_background")
    nz = len(zoom_levels)
    # The zoom strip and the cost curve get equal width, so neither dominates.
    fig = plt.figure(figsize=(6.2 + 2.15 * nz, 5.0), dpi=180)
    gs = fig.add_gridspec(1, nz + 1,
                          width_ratios=[0.80] * nz + [1.45],
                          wspace=0.20, left=0.035, right=0.985,
                          top=0.76, bottom=0.205)

    # A — the same ground at several thresholds
    zoom_axes = []
    for zi, lev in enumerate(zoom_levels):
        ax = fig.add_subplot(gs[zi])
        zoom_axes.append(ax)
        keep = zoom >= max(lev, 1e-9)
        shown = np.where(keep, zoom, np.nan)
        ax.imshow(np.ma.masked_invalid(shown), cmap="gray", vmin=0, vmax=1,
                  interpolation="nearest")
        ax.set_xticks([]); ax.set_yticks([])
        pct = 100 * keep.sum() / max((zoom > 0).sum(), 1)
        col_ = "#ff3b30" if abs(lev - PURITY_MIN) < 1e-9 else "white"
        ax.set_title(f"\u2265 {lev:g}", fontsize=11, fontweight="bold",
                     color=col_, loc="left", pad=4)
        ax.set_xlabel(f"{pct:.0f}% of cells kept", fontsize=9.5, color=col_,
                      labelpad=3)
        for sp in ax.spines.values():
            sp.set_color("#ff3b30" if col_ == "#ff3b30" else "#484f58")
        if zi == 0:
            sb = int(100 / RES)
            ax.add_patch(Rectangle((1, npx - 7), sb + 2, 6, color="black",
                                   alpha=0.75))
            ax.add_patch(Rectangle((2, npx - 3.4), sb, 1.3, color="white"))
            ax.text(2 + sb / 2, npx - 4.2, "100 m", ha="center", va="bottom",
                    fontsize=7.5, color="white")
    fig.text(0.035, 0.845, "A   The same ground, three pixel purity thresholds",
             fontsize=12, fontweight="bold", color="white", ha="left")

    # C — what the threshold costs
    ax = fig.add_subplot(gs[-1])
    col = {"intensive": "#3fb950", "extensive": "#58a6ff",
           "mown": "#d29922", "grazed": "#8250df"}
    # mown and intensive differ by under half a percent — intensive is ~95% of
    # mown by area — so a solid line simply hides the other. The broader groups
    # are drawn dashed and wider, underneath, so both stay visible.
    style = {"mown": dict(ls="--", lw=3.0, marker="", alpha=0.75, zorder=1),
             "grazed": dict(ls="--", lw=3.0, marker="", alpha=0.75, zorder=1),
             "intensive": dict(ls="-", lw=2.0, marker="o", ms=3, zorder=3),
             "extensive": dict(ls="-", lw=2.0, marker="o", ms=3, zorder=3)}
    for n in ("mown", "grazed", "intensive", "extensive"):
        if n not in curves:
            continue
        ax.plot(thresholds, curves[n], color=col[n], label=n, **style[n])
    ax.axhline(1, color="#8b949e", ls=":", lw=1.3)
    ax.axvline(PURITY_MIN, color="#ff3b30", ls="--", lw=1.8)
    ax.annotate(f"{PURITY_MIN:g}", xy=(PURITY_MIN, 0), xycoords=("data",
                "axes fraction"), xytext=(0, -17), textcoords="offset points",
                ha="center", va="top", fontsize=9.5, fontweight="bold",
                color="#ff3b30")
    ax.set_xlim(0, 0.96)
    top = max(2.0, min(4.0, max(max(v) for v in curves.values()) * 1.08))
    ax.set_ylim(0, top)
    ax.set_xlabel("purity threshold", fontsize=9.5, color="white")
    ax.set_ylabel("measured area / mapped area", fontsize=9.5, color="white")
    ax.set_title("B   What the threshold costs", fontsize=12,
                 fontweight="bold", color="white", loc="left", pad=8)
    ax.legend(fontsize=8, frameon=True, facecolor="#161b22", edgecolor="none",
              labelcolor="white", loc="upper right", ncol=1,
              bbox_to_anchor=(0.985, 0.985), framealpha=0.92,
              borderpad=0.6)

    for a in (ax,):
        a.tick_params(colors="white", labelsize=8.5)
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            a.spines[sp].set_color("#484f58")
        a.grid(alpha=0.18, ls=":", color="gray")

    fig.suptitle("Mixed pixels: why the mask keeps only field interiors",
                 fontsize=13, fontweight="bold", color="#ff3b30",
                 x=0.035, ha="left", y=0.945)
    fig.text(0.035, 0.018,
             f"Response: {psf_caption()}. "
             f"Panel A: {site} at {RES} m, {zoom_m} m across. Above 1 in "
             "panel B is the response reaching past the parcel; below 1 is "
             "edges being dropped.\nDashed = the broader groups (mown, "
             "grazed). mown and intensive nearly coincide, because intensive "
             "is most of mown by area.",
             fontsize=7.5, color="#6e7681", ha="left", va="bottom")

    out = out or (figures_dir("mask") / f"purity_explainer_{site}_{group}.png")
    fig.savefig(out, dpi=180, facecolor=fig.get_facecolor())
    print(f"  -> {out}")
    plt.show()
    return curves


def figure_mask_overview(sites: list[str] | None = None,
                         thresholds: tuple[float, ...] = (0.5, 0.9),
                         group: str = "intensive", zoom_m: float = 500,
                         rgb: dict | None = None, out: str | Path | None = None,
                         master=None, stretch: tuple[float, float] = (0.02, 0.98),
                         windows: dict | None = None, labels: dict | None = None,
                         rgb_bands: tuple[int, int, int] = (1, 2, 3),
                         verbose: bool = True) -> list | None:
    """Imagery, polygons and purity at two thresholds, one row per site.

    Reads left to right as the mask is built: what the sensor saw, what the
    cadastre says, what survives the sensor's own blur. Needs a three-band
    Sentinel-2 GeoTIFF per site at data/gee/<site>/rgb_<site>.tif; the panel says
    so if it is missing. Shows real field outlines, so real-data output is not
    publishable.

    Args:
        sites: Site keys (default: all in config.yaml).
        thresholds: Purity thresholds, one column each.
        group: Mask band.
        zoom_m: Window width in metres.
        rgb: Optional {site: path} to imagery, overriding the search.
        out: Output PNG (default: figures/mask/mask_overview_<group>.png).
        master: Master parcel set; built if None (slow).
        stretch: Percentile stretch for the imagery, per band.
        windows: Optional {site: (x0, y0, x1, y1)} in LV95, for a reproducible
            figure; otherwise the densest window is chosen.
        labels: Optional {site: row caption}.
        rgb_bands: Band indices of the imagery to show as RGB.
        verbose: Print what was found and how it was stretched.

    Returns:
        List of per-site row dicts, or None if nothing could be drawn.
    """
    import matplotlib.pyplot as plt
    import rasterio
    from matplotlib.patches import Rectangle
    from rasterio.windows import from_bounds

    from .config import CFG, figures_dir, gee_dir

    sites = list(sites or vars(CFG.sites))
    band = list(GROUPS).index(group) + 1
    npx = max(int(zoom_m / RES), 20)

    if master is None:
        print("  reading the parcel set for the polygon column...")
        master = master_parcels(verbose=False)
    codes = GROUPS[group]
    parcels = master[master["_code"].isin(codes)]

    # Locate the imagery before choosing windows. If a scene covers only part
    # of a site, the densest-parcel window can land outside it and the panel
    # comes back black — so the search is restricted to where imagery exists.
    from .config import ROOT
    scenes = {}
    for s in sites:
        p = (rgb or {}).get(s)
        if p is None:
            # Only files that announce themselves as imagery. An earlier
            # version matched any TIFF containing the site name and quietly
            # picked up a rendered z-score map instead.
            cand = (sorted(gee_dir(s).glob("rgb*.tif"))
                    + sorted((ROOT / CFG.paths.gee).glob(f"rgb*{s}*.tif"))
                    + sorted((ROOT / CFG.paths.gee).glob("rgb*.tif")))
            cand = [c for c in cand if s in c.stem or len(sites) == 1]
            p = cand[0] if cand else None
        if p is None:
            continue
        try:
            with rasterio.open(p) as d2:
                scenes[s] = (p, d2.bounds, d2.count, d2.dtypes[0])
                if verbose:
                    a = d2.read(1, out_shape=(1, min(256, d2.height),
                                              min(256, d2.width)))
                    a = a[np.isfinite(a)]
                    print(f"  {s}: {Path(p).name}, {d2.count} bands, "
                          f"{d2.dtypes[0]}, band 1 range "
                          f"{np.nanmin(a):.4g}..{np.nanmax(a):.4g}")
        except Exception as exc:
            print(f"  [!] {s}: cannot open {p} ({exc})")

    rows = []
    verify_mask()
    with rasterio.open(mask_path()) as ds:
        rb = ds.bounds
        for s in sites:
            try:
                sx0, sy0, sx1, sy1 = site_bounds_lv95(s)
            except Exception:
                print(f"  [!] unknown site {s!r}"); continue
            x0, y0 = max(rb.left, sx0), max(rb.bottom, sy0)
            x1, y1 = min(rb.right, sx1), min(rb.top, sy1)
            if s in scenes and not (windows or {}).get(s):
                sb = scenes[s][1]
                ix0, iy0 = max(x0, sb.left), max(y0, sb.bottom)
                ix1, iy1 = min(x1, sb.right), min(y1, sb.top)
                if ix1 - ix0 > zoom_m and iy1 - iy0 > zoom_m:
                    x0, y0, x1, y1 = ix0, iy0, ix1, iy1
                else:
                    print(f"  [!] {s}: the scene overlaps the site by less "
                          "than one window; using the full site extent")
            if x1 <= x0 or y1 <= y0:
                print(f"  [!] {s} lies outside the raster"); continue
            pin = (windows or {}).get(s)
            if pin:
                # An explicit window: clipped to the raster so a typo cannot
                # silently read outside it.
                px0, py0 = max(rb.left, pin[0]), max(rb.bottom, pin[1])
                px1, py1 = min(rb.right, pin[2]), min(rb.top, pin[3])
                if px1 <= px0 or py1 <= py0:
                    print(f"  [!] {s}: the given window misses the raster")
                    continue
                zoom = ds.read(band, window=from_bounds(px0, py0, px1, py1,
                                                        ds.transform))
                wx0, wy1 = px0, py1
            else:
                big = ds.read(band,
                              window=from_bounds(x0, y0, x1, y1, ds.transform))
                # densest patch, so the panel shows fields rather than empty
                # ground. Convenient, but it moves if the mask changes — pin it
                # with `windows` for a figure you need to reproduce.
                occ = (big > 0).astype(float)
                best, br, bc = -1, 0, 0
                step = max(npx // 3, 1)
                for r_ in range(0, max(big.shape[0] - npx, 1), step):
                    for c_ in range(0, max(big.shape[1] - npx, 1), step):
                        v = occ[r_:r_ + npx, c_:c_ + npx].sum()
                        if v > best:
                            best, br, bc = v, r_, c_
                zoom = big[br:br + npx, bc:bc + npx]
                wx0 = x0 + bc * RES
                wy1 = y1 - br * RES
            ext = (wx0, wx0 + zoom.shape[1] * RES,
                   wy1 - zoom.shape[0] * RES, wy1)
            rows.append({"site": s, "purity": zoom, "ext": ext})

    if not rows:
        print("nothing to draw")
        return None

    # Crisp polygon rasterisation, no PSF — the cadastre as drawn. The same
    # parcels also give the reference area for the "kept" figures: clipped to
    # the window, so the denominator is the mapped hectares actually on screen.
    # This is the definition area_check() uses, so one number means one thing
    # across the whole analysis.
    from shapely.geometry import box as _box
    cell_ha = (RES * RES) / 1e4
    for r in rows:
        e = r["ext"]
        shape = r["purity"].shape
        sub = parcels.cx[e[0]:e[1], e[2]:e[3]]
        r["poly"] = (_burn(sub.geometry, (e[0], e[3], RES), shape)
                     if len(sub) else np.zeros(shape, np.uint8))
        if len(sub):
            clip = _box(e[0], e[2], e[1], e[3])
            keep_ = sub[sub.geometry.intersects(clip)].copy()
            keep_["geometry"] = keep_.geometry.intersection(clip)
            keep_ = keep_[keep_.geometry.area > 0]
            r["parcels"] = keep_
            r["poly_ha"] = float(keep_.geometry.area.sum() / 1e4)
        else:
            r["parcels"] = sub
            r["poly_ha"] = 0.0

    # imagery, if it was exported
    for r in rows:
        r["rgb"] = None
        if r["site"] not in scenes:
            continue
        p, _b, ncount, _dt = scenes[r["site"]]
        idx = list(rgb_bands)
        if ncount < max(idx):
            print(f"  [!] {r['site']}: the file has {ncount} band(s) but "
                  f"bands {idx} were requested — export B8/B4/B3 as three "
                  "bands, or set rgb_bands")
            continue
        try:
            with rasterio.open(p) as ds:
                e = r["ext"]
                w = from_bounds(e[0], e[2], e[1], e[3], ds.transform)
                a = ds.read(indexes=idx, window=w,
                            out_shape=(3, r["purity"].shape[0],
                                       r["purity"].shape[1]),
                            boundless=True, fill_value=0).astype(np.float32)
            img = np.moveaxis(a, 0, -1)
            # Percentile stretch PER BAND on its own valid pixels. Sentinel-2
            # L2A comes as integer reflectance x10000, so a shared stretch
            # across bands would leave the image looking like a grey ramp —
            # NIR is far brighter than red or green over vegetation, and it is
            # exactly that difference the false colour is meant to show.
            for b in range(img.shape[-1]):
                ch = img[..., b]
                finite = ch[np.isfinite(ch) & (ch > 0)]
                if finite.size < 16:
                    img[..., b] = 0.0
                    continue
                lo, hi = np.quantile(finite, stretch)
                if hi - lo < 1e-9:
                    lo, hi = float(finite.min()), float(finite.max())
                img[..., b] = np.clip((ch - lo) / max(hi - lo, 1e-9), 0, 1)
            if verbose:
                mn = [float(np.nanmean(img[..., b])) for b in range(img.shape[-1])]
                print(f"  {r['site']}: stretched band means "
                      + ", ".join(f"{v:.2f}" for v in mn))
            nonzero = float((img.sum(axis=-1) > 0.01).mean())
            if nonzero < 0.5:
                print(f"  [!] {r['site']}: the exported scene covers only "
                      f"{100 * nonzero:.0f}% of the window being drawn — "
                      "re-export over the whole site, or pass a smaller "
                      "zoom_m")
            r["rgb"] = img
        except Exception as exc:
            print(f"  [!] {r['site']}: could not read {p} ({exc})")

    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Polygon as MplPoly

    plt.style.use("dark_background")
    nr, nc = len(rows), 2 + len(thresholds)
    fig, axes = plt.subplots(nr, nc, figsize=(2.9 * nc, 3.0 * nr + 0.9),
                             dpi=180, squeeze=False)
    heads = ["Sentinel-2, false colour", "Parcel polygons"] + \
            [f"Purity \u2265 {t:g}" for t in thresholds]

    for i, r in enumerate(rows):
        panels = [r["rgb"], r["poly"]] + \
                 [np.where(r["purity"] >= max(t, 1e-9), r["purity"], np.nan)
                  for t in thresholds]
        for j, dat in enumerate(panels):
            ax = axes[i][j]
            e = r["ext"]
            if j == 0:
                if dat is None:
                    ax.set_facecolor("#0d1117")
                    ax.text(0.5, 0.5, "no imagery exported\nfor this site",
                            ha="center", va="center", fontsize=8.5,
                            color="#6e7681", transform=ax.transAxes)
                else:
                    ax.imshow(dat, interpolation="nearest", extent=e)
            elif j == 1:
                # VECTOR polygons, not a raster. The point of this panel is the
                # sharpness of the boundary — rasterising it here would show the
                # same stair-stepping the purity panels have, and lose the
                # contrast the column exists to make.
                ax.set_facecolor("#000000")
                pats = []
                for geom in r["parcels"].geometry:
                    if geom is None or geom.is_empty:
                        continue
                    gs = geom.geoms if geom.geom_type == "MultiPolygon" else [geom]
                    for gg in gs:
                        pats.append(MplPoly(np.asarray(gg.exterior.coords),
                                            closed=True))
                if pats:
                    ax.add_collection(PatchCollection(
                        pats, facecolor="white", edgecolor="#9aa0a6",
                        linewidths=0.35))
                ax.set_xlim(e[0], e[1]); ax.set_ylim(e[2], e[3])
                ax.set_aspect("equal")
            else:
                ax.imshow(np.ma.masked_invalid(dat), cmap="gray", vmin=0,
                          vmax=1, interpolation="nearest", extent=e)
                keep = r["purity"] >= max(thresholds[j - 2], 1e-9)
                pct = (100 * float(keep.sum()) * cell_ha / r["poly_ha"]
                       if r["poly_ha"] > 0 else float("nan"))
                red = abs(thresholds[j - 2] - PURITY_MIN) < 1e-9
                ax.text(0.97, 0.965, f"{pct:.0f}% of mapped area",
                        transform=ax.transAxes,
                        ha="right", va="top", fontsize=9, fontweight="bold",
                        color="#ff3b30" if red else "#c9d1d9",
                        bbox=dict(facecolor="black", alpha=0.72,
                                  edgecolor="none", boxstyle="round,pad=0.25"))

            ax.set_xticks([]); ax.set_yticks([])
            if i == 0:
                ax.set_title(heads[j], fontsize=10, fontweight="bold",
                             color="white", loc="left", pad=6)
            if j == 0:
                nm = (labels or {}).get(r["site"])
                if nm is None:
                    sc = getattr(CFG.sites, r["site"], None)
                    nm = getattr(sc, "label", None) or r["site"].capitalize()
                ax.set_ylabel(nm, fontsize=11, fontweight="bold",
                              color="#ff3b30", labelpad=8)
            for sp in ax.spines.values():
                sp.set_color("#30363d")
        # scale bar on the last panel of each row
        e = r["ext"]
        w_ = e[1] - e[0]
        bx, by = e[0] + 0.05 * w_, e[2] + 0.05 * (e[3] - e[2])
        axL = axes[i][-1]
        axL.add_patch(Rectangle((bx - 0.012 * w_, by - 0.012 * w_),
                                100 + 0.024 * w_, 0.062 * w_,
                                color="black", alpha=0.72, zorder=5))
        axL.add_patch(Rectangle((bx, by), 100, 0.009 * w_, color="white",
                                zorder=6))
        axL.text(bx + 50, by + 0.016 * w_, "100 m", ha="center", va="bottom",
                 fontsize=7, color="white", zorder=6)

    fig.suptitle("From imagery to mask: what each step keeps",
                 fontsize=12.5, fontweight="bold", color="#ff3b30",
                 x=0.012, ha="left", y=0.985)
    fig.supxlabel(f"{group} grassland at {RES} m, {zoom_m} m across. Purity is "
                  "the share of a cell's sensor response coming from inside a "
                  f"parcel ({psf_caption()}). White "
                  "is pure, grey is mixed, black is discarded.\n"
                  "Percentages are retained area against the polygon area in "
                  "the same window \u2014 above 100% at a low threshold means "
                  "the response is reaching past the parcels.",
                  fontsize=7.5, color="#6e7681")
    fig.tight_layout(rect=[0.004, 0.055, 0.998, 0.958])
    fig.subplots_adjust(hspace=0.06, wspace=0.035)
    out = out or (figures_dir("mask") / f"mask_overview_{group}.png")
    fig.savefig(out, dpi=180, facecolor=fig.get_facecolor())
    print(f"  -> {out}")
    plt.show()
    return rows


def figure_parcel_size(master=None, group: str | None = None, out: str | Path | None = None,
                       year: int | None = None,
                       bins: tuple[float, ...] = (0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0),
                       shape_clip: float = 16) -> dict | None:
    """How big and how convoluted the fields are: area and shape index per category.

    The shape index (perimeter / sqrt(area), 4.0 = square) explains what the mask
    costs better than area alone: a long, narrow field is mostly edge however
    many pixels fit inside it. Distributions only, so publishable for real data.

    Args:
        master: Master parcel set; built if None (slow).
        group: One category, or None for all four.
        out: Output PNG (default: figures/mask/parcel_size.png).
        year: Unused; kept for call compatibility.
        bins: Area tick marks on the log axis, in ha.
        shape_clip: Shape index above which the right panel clips, for display only.

    Returns:
        dict category -> (areas in ha, shape indices), or None if no category has
        at least 10 parcels.
    """
    import matplotlib.pyplot as plt

    from .config import figures_dir

    if master is None:
        print("  building the parcel set — pass master= to skip this")
        master = master_parcels(verbose=False)
    groups = [group] if group else ["intensive", "extensive", "mown", "grazed"]

    data = {}
    for n in groups:
        sub = master[master["_code"].isin(GROUPS[n])]
        if len(sub) < 10:
            continue
        a = (sub.geometry.area / 1e4).to_numpy()
        shp = (sub.geometry.length / np.sqrt(sub.geometry.area)).to_numpy()
        ok = np.isfinite(a) & np.isfinite(shp) & (a > 0)
        data[n] = (a[ok], shp[ok])
    if not data:
        print("no group has enough parcels")
        return None
    names = list(data)
    col = {"intensive": "#3fb950", "extensive": "#58a6ff",
           "mown": "#d29922", "grazed": "#8250df"}

    plt.style.use("dark_background")
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.2), dpi=180,
                             constrained_layout=True)

    ax = axes[0]
    vals = [np.log10(data[n][0]) for n in names]
    parts = ax.violinplot(vals, positions=range(len(names)), widths=0.85,
                          showextrema=False, showmedians=True)
    for b, n, v in zip(parts["bodies"], names, vals):
        b.set_facecolor(col[n]); b.set_alpha(0.78); b.set_edgecolor("none")
        vt = b.get_paths()[0].vertices
        vt[:, 1] = np.clip(vt[:, 1], v.min(), v.max())
    parts["cmedians"].set_color("white"); parts["cmedians"].set_linewidth(1.9)
    ax.set_yticks([np.log10(b) for b in bins])
    ax.set_yticklabels([f"{b:g}" for b in bins])
    ax.set_ylabel("parcel area [ha], log scale", fontsize=10, color="white")
    ax.set_title("How big is a field?", fontsize=11.5, fontweight="bold",
                 color="white", loc="left")
    # what a 10 m pixel is, on the same axis
    ax.axhline(np.log10(0.01), color="#eda100", ls="--", lw=1.5)
    ax.annotate("one 10 m pixel = 0.01 ha", xy=(0.985, np.log10(0.01)),
                xycoords=("axes fraction", "data"), xytext=(0, 6),
                textcoords="offset points", ha="right", fontsize=9,
                color="#eda100")

    ax = axes[1]
    # The shape index has a long upper tail — a minority of parcels run well
    # past 10 — so the violins are clipped for display. Clipping too low makes
    # every violin end in a flat top, which reads as data at the ceiling rather
    # than as a cut-off; the axis says so when the clip actually bites.
    n_clipped = int(sum(int((data[n][1] > shape_clip).sum()) for n in names))
    vals2 = [np.clip(data[n][1], 3.6, shape_clip) for n in names]
    parts = ax.violinplot(vals2, positions=range(len(names)), widths=0.85,
                          showextrema=False, showmedians=True)
    for b, n, v in zip(parts["bodies"], names, vals2):
        b.set_facecolor(col[n]); b.set_alpha(0.78); b.set_edgecolor("none")
        vt = b.get_paths()[0].vertices
        vt[:, 1] = np.clip(vt[:, 1], v.min(), v.max())
    parts["cmedians"].set_color("white"); parts["cmedians"].set_linewidth(1.9)
    ax.axhline(4.0, color="#eda100", ls="--", lw=1.5)
    ax.annotate("4.0 = a perfect square", xy=(0.985, 4.0),
                xycoords=("axes fraction", "data"), xytext=(0, 6),
                textcoords="offset points", ha="right", fontsize=9,
                color="#eda100")
    ax.set_ylabel("perimeter / \u221aarea" +
                  (f"  (display clipped at {shape_clip:g})" if n_clipped else ""),
                  fontsize=10, color="white")
    ax.set_title("How much of it is edge?", fontsize=11.5, fontweight="bold",
                 color="white", loc="left")

    for ax_ in axes:
        ax_.set_xticks(range(len(names)))
        ax_.set_xticklabels(
            [f"{n}\nn={len(data[n][0]):,}" for n in names], fontsize=9)
        ax_.grid(alpha=0.22, ls=":", color="gray")
        ax_.tick_params(colors="white", labelsize=9)
        for sp in ("top", "right"):
            ax_.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax_.spines[sp].set_color("#484f58")

    print(f"  {'group':<12}{'parcels':>9}{'median ha':>11}{'<0.5 ha':>9}"
          f"{'median shape':>14}{'pixels in a median field':>26}")
    for n in names:
        a, shp = data[n]
        print(f"  {n:<12}{len(a):>9,}{np.median(a):>11.2f}"
              f"{100 * (a < 0.5).mean():>8.0f}%{np.median(shp):>14.2f}"
              f"{np.median(a) / 0.01:>26.0f}")
    if n_clipped:
        tot = sum(len(data[n][0]) for n in names)
        print(f"\n  {n_clipped:,} of {tot:,} parcels ({100 * n_clipped / tot:.1f}%) "
              f"have a shape index above {shape_clip:g} and are clipped in the "
              "right panel only.")

    fig.suptitle("Swiss grassland parcels: small, and mostly edge",
                 fontsize=12.5, fontweight="bold", color="#ff3b30",
                 x=0.012, ha="left")
    # Measured, not written in: the caption once claimed "about a hundred
    # pixels" long after the size filter that made it true had gone.
    cells = " · ".join(f"{n} {np.median(data[n][0]) / 0.01:.0f}" for n in names)
    fig.supxlabel(f"10 m cells in a median field:  {cells}", fontsize=7.5,
                  color="#6e7681")
    out = out or (figures_dir("mask") / "parcel_size.png")
    fig.savefig(out, dpi=180, facecolor=fig.get_facecolor())
    print(f"\n  -> {out}")
    plt.show()
    return data


def figure_psf(out: str | Path | None = None, mtf: float | None = None,
               span: float = 40, ground: bool = True) -> dict:
    """What one pixel actually sees: the modelled response, in section and as a surface.

    The detector box convolved with Gaussian optics, derived from the system MTF
    (METHOD §2.3). Uses no data, so it is rendered unstamped by `make gallery`.

    Args:
        out: Output PNG (default: figures/mask/psf.png).
        mtf: System MTF at Nyquist (default: config).
        span: Metres of ground shown.
        ground: Draw the ground plane under the surface.

    Returns:
        dict with sigma, fwhm, mtf_detector, mtf_optics.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as MplPoly
    from scipy.special import erf

    from .config import figures_dir

    pp = psf_params(mtf=mtf, model="box_gauss")
    m_, mtf_det, mtf_opt = pp["mtf_system"], pp["mtf_detector"], pp["mtf_optics"]
    sig, fwhm = pp["sigma_m"], pp["fwhm_m"]

    def prof(x):
        return 0.5 * (erf((x + RES / 2) / (sig * np.sqrt(2)))
                      - erf((x - RES / 2) / (sig * np.sqrt(2))))

    pk = prof(0.0)
    xs = np.linspace(-span / 2, span / 2, 1200)
    ys = prof(xs) / pk

    plt.style.use("dark_background")
    fig = plt.figure(figsize=(13.4, 6.4), dpi=180)
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.18], wspace=0.10,
                          left=0.055, right=0.985, top=0.87, bottom=0.20)

    # ---- cross-section ---------------------------------------------------
    ax = fig.add_subplot(gs[0])
    ax.axvspan(-RES / 2, RES / 2, color="#eda100", alpha=0.13, lw=0)
    for e in (-RES / 2, RES / 2):
        ax.axvline(e, color="#eda100", ls=":", lw=1.2)
    ax.fill_between(xs, ys, color="white", alpha=0.10, lw=0)
    ax.plot(xs, ys, color="white", lw=2.6)
    ax.plot([-fwhm / 2, fwhm / 2], [0.5, 0.5], color="#eda100", lw=1.8)
    for e in (-fwhm / 2, fwhm / 2):
        ax.plot([e, e], [0.46, 0.54], color="#eda100", lw=1.8)
    ax.annotate(f"FWHM {fwhm:.1f} m", xy=(fwhm / 2, 0.5),
                xytext=(fwhm / 2 + 2.5, 0.60), fontsize=10, color="#eda100",
                arrowprops=dict(arrowstyle="-", color="#eda100", lw=0.9))
    ax.annotate(f"the nominal {RES:g} m pixel", xy=(RES / 2, 0.93),
                xytext=(RES / 2 + 3, 1.02), fontsize=10, color="#eda100",
                arrowprops=dict(arrowstyle="-", color="#eda100", lw=0.9))
    ax.set_xlim(-span / 2, span / 2); ax.set_ylim(0, 1.16)
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_xlabel("metres from the pixel centre", fontsize=10, color="white")
    ax.set_ylabel("response, normalised to peak", fontsize=10, color="white")
    ax.set_title("In cross-section", fontsize=11.5, fontweight="bold",
                 color="white", loc="left")
    ax.grid(alpha=0.18, ls=":", color="gray")

    # ---- surface over ground --------------------------------------------
    ax2 = fig.add_subplot(gs[1])
    ax2.set_axis_off()
    yaw, pit, S = -0.70, 0.56, 1.0
    cy, sy = np.cos(yaw), np.sin(yaw)
    cp, sp = np.cos(pit), np.sin(pit)

    def P(u, v, w):
        # NOTE the sign on the vertical: matplotlib's y axis points up while a
        # canvas points down. Porting the projection without flipping it put
        # the surface in a pit and the ground plane in the sky.
        x, y, z = u - 0.5, v - 0.5, w - 0.5
        return ((x * cy - y * sy) * 0.82 * S,
                -((x * sy + y * cy) * cp - z * sp) * 0.42 * S)

    def draw(pts, fc=None, ec=None, lw=1.0, z=1):
        ax2.add_patch(MplPoly(np.array(pts), closed=True,
                              facecolor=fc or "none", edgecolor=ec or "none",
                              lw=lw, zorder=z))

    GZ = -1.95 if ground else 0.0
    if ground:
        gq = [P(0, 0, GZ), P(1, 0, GZ), P(1, 1, GZ), P(0, 1, GZ)]
        draw(gq, fc="#3d4a32", z=1)
        bnd = [0.0, 0.26, 0.50, 0.74, 1.0]
        tint = ["#6f8348", "#93a05c", "#b0a068", "#5c7340"]
        for k in range(4):
            draw([P(0, bnd[k], GZ), P(1, bnd[k] - 0.08, GZ),
                  P(1, bnd[k + 1] - 0.08, GZ), P(0, bnd[k + 1], GZ)],
                 fc=tint[k], z=2)
        uu = np.linspace(0, 1, 12)
        rv = 0.66 - 0.26 * uu
        riv = [P(u, r + 0.028, GZ) for u, r in zip(uu, rv)] + \
              [P(u, r - 0.028, GZ) for u, r in zip(uu[::-1], rv[::-1])]
        draw(riv, fc="#46687c", z=3)
        draw(gq, ec="#ffffff55", lw=1.2, z=4)

    hp = (RES / span) / 2
    sq = [(0.5 - hp, 0.5 - hp), (0.5 + hp, 0.5 - hp),
          (0.5 + hp, 0.5 + hp), (0.5 - hp, 0.5 + hp)]
    if ground:
        draw([P(u, v, GZ) for u, v in sq], fc="#eda10033", ec="#eda100",
             lw=1.8, z=5)
        for u, v in sq:
            a, b = P(u, v, 0), P(u, v, GZ)
            ax2.plot([a[0], b[0]], [a[1], b[1]], color="#eda100", lw=0.9,
                     ls=(0, (3, 3)), zorder=5)

    draw([P(0, 0, 0), P(1, 0, 0), P(1, 1, 0), P(0, 1, 0)],
         fc="#ffffff08", ec="#ffffff6b", lw=1.3, z=6)
    for t in (0.25, 0.5, 0.75):
        for a, b in ((P(t, 0, 0), P(t, 1, 0)), (P(0, t, 0), P(1, t, 0))):
            ax2.plot([a[0], b[0]], [a[1], b[1]], color="#ffffff33", lw=0.9,
                     zorder=6)

    N = 30
    g1 = np.linspace(-span / 2, span / 2, N + 1)
    f1 = prof(g1) / pk
    quads = []
    for i in range(N):
        for j in range(N):
            pts = [P(i / N, j / N, f1[i] * f1[j]),
                   P((i + 1) / N, j / N, f1[i + 1] * f1[j]),
                   P((i + 1) / N, (j + 1) / N, f1[i + 1] * f1[j + 1]),
                   P(i / N, (j + 1) / N, f1[i] * f1[j + 1])]
            d = ((i / N - 0.5) * sy + (j / N - 0.5) * cy) * sp \
                + (f1[i] * f1[j] - 0.5) * cp
            quads.append((d, pts))
    for _d, pts in sorted(quads, key=lambda q: q[0]):
        draw(pts, fc="#000000", ec="#ffffffd0", lw=0.45, z=7)

    draw([P(u, v, 0) for u, v in sq], fc="#eda10055", ec="#eda100", lw=2,
         z=8)
    c = P(0.5 + hp, 0.5 - hp, 0)
    ax2.annotate(f"{RES:g} \u00d7 {RES:g} m pixel", xy=c,
                 xytext=(c[0] + 0.12, c[1] - 0.10), fontsize=9.5,
                 color="#eda100",
                 arrowprops=dict(arrowstyle="-", color="#eda100", lw=0.9))
    # Edge labels offset outward from the floor centre, so they never land on
    # the mesh.
    fc = P(0.5, 0.5, 0)
    for (a_, b_) in ((P(0, 0, 0), P(1, 0, 0)), (P(1, 0, 0), P(1, 1, 0))):
        mx, my = (a_[0] + b_[0]) / 2, (a_[1] + b_[1]) / 2
        dx, dy = mx - fc[0], my - fc[1]
        n = max((dx * dx + dy * dy) ** 0.5, 1e-6)
        ax2.plot([a_[0] + dx / n * 0.05, b_[0] + dx / n * 0.05],
                 [a_[1] + dy / n * 0.05, b_[1] + dy / n * 0.05],
                 color="#ffffff70", lw=1, zorder=9)
        ax2.text(mx + dx / n * 0.12, my + dy / n * 0.12, f"{span:g} m",
                 fontsize=9, color="#ffffffb0", ha="center", va="center",
                 zorder=9)
    if ground:
        gl = P(0, 0.5, GZ)
        ax2.text(gl[0] - 0.03, gl[1], "the ground\nbelow", fontsize=9.5,
                 color="#ffffffcc", ha="right", va="center")
    ax2.set_title(f"Over {span:g} \u00d7 {span:g} m of ground",
                  fontsize=11.5, fontweight="bold", color="white", loc="left",
                  pad=12)
    # Limits fitted to the content, so the panel is not shrunk by empty space.
    ax2.set_xlim(-0.66, 0.78)
    ax2.set_ylim(-0.86 if ground else -0.30, 0.30)
    ax2.set_aspect("equal")

    for a_ in (ax,):
        a_.tick_params(colors="white", labelsize=9)
        for sp_ in ("top", "right"):
            a_.spines[sp_].set_visible(False)
        for sp_ in ("left", "bottom"):
            a_.spines[sp_].set_color("#484f58")

    vals = [("\u03c3 optics", f"{sig:.2f} m"), ("FWHM", f"{fwhm:.1f} m"),
            ("detector MTF", f"{mtf_det:.3f}"),
            ("optics MTF", f"{mtf_opt:.3f}"), ("system MTF", f"{m_:.2f}")]
    for i, (k, v) in enumerate(vals):
        x = 0.055 + i * 0.185
        fig.text(x, 0.115, k, fontsize=8.5, color="#8b949e")
        fig.text(x, 0.072, v, fontsize=15, fontweight="bold", color="white")

    fig.suptitle("A pixel is blurrier than its own grid", fontsize=13,
                 fontweight="bold", color="#ff3b30", x=0.012, ha="left",
                 y=0.965)
    fig.text(0.055, 0.018,
             "MTF at Nyquist. The detector footprint, a box one pixel wide, "
             "convolved with a Gaussian optical response \u2014 which is why "
             "the top is flatter than a Gaussian.\nApproximations: the optics "
             "are Airy-like rather than Gaussian, the real response is "
             "slightly wider along-track, and one representative MTF stands "
             "for all bands.", fontsize=7, color="#6e7681")
    out = out or (figures_dir("mask") / "psf.png")
    fig.savefig(out, dpi=180, facecolor=fig.get_facecolor())
    print(f"  sigma {sig:.2f} m, FWHM {fwhm:.1f} m  ->  {out}")
    plt.show()
    return {"sigma": float(sig), "fwhm": float(fwhm),
            "mtf_detector": float(mtf_det), "mtf_optics": float(mtf_opt)}
