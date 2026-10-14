"""Figures.

Output goes to output/figures/<site>/. Every figure names its site, window and
regime, so a screenshot is readable without the notebook beside it. Dark ground,
grey baseline, white dashed median, red current year — the same grammar as the
precipitation panels in the deck.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .config import CFG, figures_dir
from .detect import detect_all, detect_parcel_year
from .metrics import equalised, last_cut, last_cut_paired, late_fraction

REGIME_TITLE = {
    "intensive": "Intensively used permanent meadows  ·  LNF 613",
    "extensive": "Extensive & wenig intensiv meadows  ·  LNF 611/612",
}
REGIME_NOTE = {
    "intensive": "4-6 cuts per season, first cut early May, no legal "
                 "restriction on the cutting date",
    "extensive": "~2 cuts per season, cutting not permitted before 15 June "
                 "in the valley and hill zones",
}


def _style(axes, square=False):
    for a in np.atleast_1d(axes).ravel():
        if square:
            a.set_box_aspect(1)
        a.grid(alpha=0.22, ls=":", color="gray")
        a.tick_params(colors="white", labelsize=9, labelbottom=True)
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            a.spines[s].set_color("#484f58")


def _save(fig, site, name):
    p = figures_dir(site) / f"{name}.png"
    fig.savefig(p, dpi=180, facecolor=fig.get_facecolor())
    print(f"  -> {p}")
    return p


def figure_last_cut(df, site, regime="intensive", cfg=None):
    """THE headline figure: when cutting stopped.

    The most robust of the cut metrics — it needs the detector to find one cut
    per parcel-year rather than all of them, and the observational bias runs the
    right way: the drought year has the most images, so a late cut is easier to
    detect there, not harder.
    """
    C = cfg or CFG
    cur = C.project.current_year
    lc = last_cut(df, regime, C)
    if lc.empty:
        print("no parcel-years with a detected cut")
        return None
    yrs = sorted(lc.year.unique())
    base = [y for y in yrs if y != cur]

    plt.style.use("dark_background")
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.2), dpi=180,
                             constrained_layout=True)

    ax = axes[0]
    data = [lc.loc[lc.year == y, "last_cut"].to_numpy() for y in yrs]
    bp = ax.boxplot(data, positions=range(len(yrs)), widths=0.6,
                    patch_artist=True, medianprops=dict(color="white", lw=2),
                    flierprops=dict(marker=".", markersize=3,
                                    markerfacecolor="#8b949e",
                                    markeredgecolor="none"))
    for patch, y in zip(bp["boxes"], yrs):
        patch.set_facecolor("#ff3b30" if y == cur else "#8b949e")
        patch.set_alpha(0.85); patch.set_edgecolor("none")
    for w in bp["whiskers"] + bp["caps"]:
        w.set_color("#8b949e")
    bmed = lc[lc.year.isin(base)].last_cut.median()
    ax.axhline(bmed, color="white", ls="--", lw=1.4,
               label=f"baseline median DOY {bmed:.0f}")
    ax.set_xticks(range(len(yrs))); ax.set_xticklabels(yrs)
    ax.set_ylabel("day of year of the LAST detected cut", fontsize=10,
                  color="white")
    ax.set_title("When cutting stopped", fontsize=11.5, fontweight="bold",
                 color="white", loc="left")
    ax.legend(fontsize=8, frameon=True, facecolor="#161b22", edgecolor="none",
              labelcolor="white", loc="lower left")

    ax = axes[1]
    pair = last_cut_paired(lc, C)
    jit = np.random.default_rng(0).normal(0, 1.2, (len(pair), 2))
    ax.scatter(pair.baseline + jit[:, 0], pair.current + jit[:, 1], s=18,
               color="#ff3b30", alpha=0.45, edgecolor="none")
    lo = float(min(pair.baseline.min(), pair.current.min())) - 5
    hi = float(max(pair.baseline.max(), pair.current.max())) + 5
    ax.plot([lo, hi], [lo, hi], color="white", ls="--", lw=1.4,
            label="no change")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel("parcel's own baseline mean", fontsize=10, color="white")
    ax.set_ylabel(f"same parcel in {cur}", fontsize=10, color="white")
    ax.set_title("Each parcel against itself", fontsize=11.5,
                 fontweight="bold", color="white", loc="left")
    ax.text(0.03, 0.97,
            f"median shift {pair['shift_days'].median():+.0f} days\n"
            f"{100 * (pair['shift_days'] < 0).mean():.0f}% of {len(pair)} parcels "
            "stopped earlier",
            transform=ax.transAxes, va="top", fontsize=9, color="#c9d1d9",
            bbox=dict(facecolor="#161b22", edgecolor="#30363d",
                      boxstyle="round,pad=0.4"))
    ax.legend(fontsize=8, frameon=True, facecolor="#161b22", edgecolor="none",
              labelcolor="white", loc="lower right")

    _style(axes, square=True)
    fig.suptitle(f"{REGIME_TITLE[regime]}   ·   {site}   ·   last cut",
                 fontsize=11.5, fontweight="bold", color="#ff3b30")
    fig.supxlabel("All years truncated at the current year's last observation. "
                  "Needs only ONE cut found per parcel-year, so far less "
                  "sensitive to detection rate than a count or a share.",
                  fontsize=7, color="#6e7681")
    _save(fig, site, f"last_cut_{regime}")
    plt.show()
    return lc


def figure_summary(df, site, regime="intensive", cfg=None, n_rep=10):
    """Late-season share, with the reason counts cannot be used beside it.

    The bottom panel is NOT a result. Raw counts follow image availability;
    matched to an equal image count the years are indistinguishable. Showing
    both is what justifies using a ratio in the panel above.
    """
    C = cfg or CFG
    cur = C.project.current_year
    split = C.analysis.late_split_doy
    overall, _, _ = late_fraction(df, regime, C, split)
    ev = detect_all(df, C, verbose=False)
    codes = C.landuse.intensive if regime == "intensive" else C.landuse.extensive
    ev = ev[ev.lnf_code.isin(codes)]
    counts = ev.groupby("year").agg(cuts=("n_cuts", "mean"),
                                    obs=("n_obs", "mean"))
    eq = equalised(df, regime, C, n_rep=n_rep)["mean_cuts"]

    yrs = list(overall.index)
    base = [y for y in yrs if y != cur]
    col = ["#ff3b30" if y == cur else "#8b949e" for y in yrs]
    x = np.arange(len(yrs))

    plt.style.use("dark_background")
    fig, axes = plt.subplots(2, 1, figsize=(7.4, 7.8), dpi=180, sharex=True,
                             constrained_layout=True)

    ax = axes[0]
    ax.axhspan(overall.loc[base, "frac_late"].min(),
               overall.loc[base, "frac_late"].max(), color="#8b949e",
               alpha=0.14, lw=0, zorder=0, label="baseline range")
    ax.bar(x, overall.frac_late.reindex(yrs), color=col, width=0.68, zorder=2)
    bmed = overall.loc[base, "frac_late"].median()
    ax.axhline(bmed, color="white", ls="--", lw=1.6, zorder=3,
               label=f"baseline median {bmed:.0%}")
    for xi, y in enumerate(yrs):
        v = overall.frac_late[y]
        ax.text(xi, v + 0.015, f"{v:.0%}", ha="center", fontsize=8.5,
                color="#ff3b30" if y == cur else "#c9d1d9")
    ax.set_ylim(0, max(overall.frac_late) * 1.30)
    ax.set_ylabel(f"share of cuts after DOY {split}", fontsize=10, color="white")
    ax.set_title("Late-season mowing", fontsize=11.5, fontweight="bold",
                 color="white", loc="left")
    ax.legend(fontsize=7.5, frameon=True, facecolor="#161b22",
              edgecolor="none", labelcolor="white", loc="upper left", ncol=2)

    ax = axes[1]
    w = 0.34
    ax.bar(x - w / 2, counts.cuts.reindex(yrs), width=w, color=col, alpha=0.35,
           edgecolor="#8b949e", lw=0.8, label="raw (confounded)")
    ax.bar(x + w / 2, eq.reindex(yrs), width=w, color=col,
           label="matched to equal image count")
    ax2 = ax.twinx()
    ax2.plot(x, counts.obs.reindex(yrs), "o-", color="#58a6ff", lw=1.8, ms=5,
             label="clear observations")
    ax2.set_ylabel("clear observations per parcel", fontsize=9, color="#58a6ff")
    ax2.tick_params(colors="#58a6ff", labelsize=8)
    ax2.spines["top"].set_visible(False)
    ax.set_ylabel("mean cuts per parcel", fontsize=10, color="white")
    ax.set_xticks(x); ax.set_xticklabels(yrs)
    ax.set_ylim(0, max(counts.cuts.max(), eq.max()) * 1.35)
    ax.set_title("Cut counts: raw vs. matched observing conditions",
                 fontsize=11.5, fontweight="bold", color="white", loc="left")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=7.5, frameon=True,
              facecolor="#161b22", edgecolor="none", labelcolor="white",
              loc="upper left", ncol=3)

    _style(axes)
    n_parcels = df[df.lnf_code.isin(codes)].parcel_id.nunique()
    fig.suptitle(f"{REGIME_TITLE[regime]}   ·   {site}   "
                 f"(n = {n_parcels} parcels)",
                 fontsize=11.5, fontweight="bold", color="#ff3b30")
    fig.supxlabel(f"{REGIME_NOTE[regime]}\nSentinel-2 {C.detect.index}  ·  "
                  "shares, not counts: detection depends on how many clear "
                  "images a year had", fontsize=7, color="#6e7681")
    _save(fig, site, f"mowing_summary_{regime}")
    plt.show()
    return overall.join(counts).join(eq.rename("cuts_matched"))


def figure_legal_check(df, site, cfg=None):
    """Validation, not a result.

    Extensive meadows may not be cut before 15 June in the valley and hill
    zones, so detections left of that line are known errors. This turns a
    regulation into ground truth — no synthetic data, no assumptions.
    """
    C = cfg or CFG
    legal = C.landuse.earliest_legal_cut_doy
    cur = C.project.current_year
    de = df[df.lnf_code.isin(C.landuse.extensive)]
    if de.empty:
        print("no extensive parcels in this sample")
        return None
    rows = []
    for (_pid, yr), g in de.groupby(["parcel_id", "year"]):
        for d in (detect_parcel_year(g, C)["cut_dates"] or []):
            rows.append({"year": yr, "doy": d, "cur": yr == cur})
    cuts = pd.DataFrame(rows)
    if cuts.empty:
        print("no cuts detected on extensive parcels")
        return None

    plt.style.use("dark_background")
    fig, ax = plt.subplots(figsize=(7.6, 4.8), dpi=180, constrained_layout=True)
    ax.scatter(cuts.doy[~cuts.cur], cuts.year[~cuts.cur], s=26,
               color="#8b949e", alpha=0.75, label="baseline")
    ax.scatter(cuts.doy[cuts.cur], cuts.year[cuts.cur], s=34,
               color="#ff3b30", label=str(cur))
    ax.axvline(legal, color="#58a6ff", lw=1.8)
    ax.axvspan(C.detect.season[0], legal, color="#58a6ff", alpha=0.10, lw=0)
    ax.annotate("15 June\ncutting not permitted before",
                xy=(legal - 3, 0.02), xycoords=("data", "axes fraction"),
                ha="right", va="bottom", fontsize=8, color="#58a6ff")
    bad = int((cuts.doy < legal).sum())
    ax.set_title(f"Detector check on extensive meadows — {site}",
                 fontsize=11, fontweight="bold", color="white", loc="left")
    ax.text(0.0, 1.02, f"{bad} of {len(cuts)} detections "
            f"({100 * bad / len(cuts):.0f}%) fall in the legally prohibited "
            "period", transform=ax.transAxes, fontsize=8.5, color="#c9d1d9",
            va="bottom")
    ax.set_xlabel("day of year", fontsize=10, color="white")
    ax.set_ylabel("year", fontsize=10, color="white")
    ax.set_xlim(*C.detect.season)
    ax.legend(fontsize=8, frameon=True, facecolor="#161b22", edgecolor="none",
              labelcolor="white", loc="upper right")
    _style(ax)
    _save(fig, site, "mowing_legal_check")
    plt.show()
    return cuts


def plot_parcel(df, pid, ev, site=None, years=None, cfg=None, save=False):
    """One parcel across years with detected cuts marked.

    Look at a handful before believing any count. The question for each panel:
    is there a red line at every place your eye sees an abrupt drop, and
    nowhere else?
    """
    C = cfg or CFG
    years = years if years is not None else sorted(
        df.loc[df.parcel_id == pid, "year"].unique())
    fig, axes = plt.subplots(len(years), 1, figsize=(10, 1.9 * len(years)),
                             sharex=True, sharey=True, constrained_layout=True)
    axes = np.atleast_1d(axes)
    plt.style.use("default")
    for ax, y in zip(axes, years):
        g = df[(df.parcel_id == pid) & (df.year == y)].sort_values("doy")
        ok = g[~g.haze] if "haze" in g else g
        ax.plot(ok.doy, ok[C.detect.index], "o-", ms=3.5, lw=1.1,
                color="#2166ac", label=C.detect.index, zorder=3)
        ax.plot(ok.doy, ok[C.detect.confirm], "o-", ms=3, lw=1.0,
                color="#b2182b", label=C.detect.confirm, zorder=3)
        if "haze" in g and g.haze.any():
            hz = g[g.haze]
            ax.plot(hz.doy, hz[C.detect.index], "x", ms=6, color="#999999",
                    label="flagged haze", zorder=4)
        r = ev[(ev.parcel_id == pid) & (ev.year == y)]
        title = str(y)
        if len(r):
            row = r.iloc[0]
            for c in (row.cut_dates or []):
                ax.axvline(c, color="#ff3b30", lw=1.4, alpha=0.75, zorder=2)
            n = row.n_cuts
            title = (f"{y}   {'' if pd.isna(n) else int(n)} cuts   "
                     f"threshold {row.thr:.2f}   max gap {row.max_gap:.0f} d")
        ax.set_title(title, fontsize=8.5, loc="left")
        ax.grid(alpha=0.3, ls=":")
        ax.set_ylim(0, 1)
    axes[0].legend(fontsize=7.5, loc="lower left", ncol=3)
    axes[-1].set_xlabel("day of year")
    fig.suptitle(f"parcel {pid}", fontsize=10)
    if save and site:
        _save(fig, site, f"parcel_{str(pid)[-8:]}")
    plt.show()


def figure_cut_timing(df, site, regime="intensive", cfg=None, show_split=False):
    """Distribution of cut dates per year — no arbitrary threshold.

    Replaces the DOY-190 bar chart. That split had to be chosen, and any choice
    invites the question "why 190?". Showing the whole distribution answers it
    by not asking: the reader sees where the cuts actually fell.

    Left  — violin per year, all detected cuts pooled across parcels. The SHAPE
            is the message; a violin is normalised, so it is not distorted by
            one year having more detections than another.
    Right — cumulative share of each year's cuts by day. A line that rises
            earlier means cutting was front-loaded. This is the panel to read
            for a shift in timing: it is a ratio at every point, so it survives
            the fact that detection rate differs between years.
    """
    C = cfg or CFG
    cur = C.project.current_year
    codes = C.landuse.intensive if regime == "intensive" else C.landuse.extensive
    d = df[df.lnf_code.isin(codes)] if "lnf_code" in df else df
    if C.detect.truncate_to_current:
        from .detect import season_end
        d = d[d.doy <= season_end(d, C)]

    rows = []
    for (pid, yr), g in d.groupby(["parcel_id", "year"]):
        for c in (detect_parcel_year(g, C)["cut_dates"] or []):
            rows.append({"year": int(yr), "doy": float(c), "parcel_id": pid})
    cuts = pd.DataFrame(rows)
    if cuts.empty:
        print("no cuts detected")
        return None
    yrs = sorted(cuts.year.unique())
    base = [y for y in yrs if y != cur]

    plt.style.use("dark_background")
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.4), dpi=180,
                             constrained_layout=True)

    ax = axes[0]
    data = [cuts.loc[cuts.year == y, "doy"].to_numpy() for y in yrs]
    parts = ax.violinplot(data, positions=range(len(yrs)), widths=0.85,
                          showextrema=False, showmedians=True)
    for b, y, v in zip(parts["bodies"], yrs, data):
        b.set_facecolor("#ff3b30" if y == cur else "#8b949e")
        b.set_alpha(0.85 if y == cur else 0.55)
        b.set_edgecolor("none")
        # CLIP the kernel tails to the observed range. A violin's KDE extends
        # past the data, and here that draws density beyond the truncation day
        # — where a cut cannot exist, because there are no images. Cutting the
        # polygon at the data range stops the figure implying otherwise.
        verts = b.get_paths()[0].vertices
        verts[:, 1] = np.clip(verts[:, 1], v.min(), v.max())
    parts["cmedians"].set_color("white")
    parts["cmedians"].set_linewidth(1.8)
    # the individual cuts, so the reader sees the sample behind the shape
    for i, v in enumerate(data):
        jit = np.random.default_rng(i).normal(0, 0.055, len(v))
        ax.plot(np.full(len(v), i) + jit, v, ".", ms=1.8,
                color="#ffffff", alpha=0.25, zorder=3)
    ax.set_xticks(range(len(yrs))); ax.set_xticklabels(yrs)
    ax.set_ylabel("day of year of each detected cut", fontsize=10,
                  color="white")
    ax.set_title("When cuts happened", fontsize=11.5, fontweight="bold",
                 color="white", loc="left")
    if show_split:
        ax.axhline(C.analysis.late_split_doy, color="#58a6ff", ls=":", lw=1.2)

    ax = axes[1]
    grid = np.arange(C.detect.season[0], C.detect.season[1] + 1)
    for y in yrs:
        v = np.sort(cuts.loc[cuts.year == y, "doy"].to_numpy())
        ecdf = np.searchsorted(v, grid, side="right") / len(v)
        ax.plot(grid, ecdf, lw=2.4 if y == cur else 1.0,
                color="#ff3b30" if y == cur else "#8b949e",
                alpha=1.0 if y == cur else 0.55,
                label=str(y) if y == cur else None, zorder=3 if y == cur else 1)
    base_stack = []
    for y in base:
        v = np.sort(cuts.loc[cuts.year == y, "doy"].to_numpy())
        base_stack.append(np.searchsorted(v, grid, side="right") / len(v))
    bs = np.vstack(base_stack)
    ax.fill_between(grid, bs.min(axis=0), bs.max(axis=0), color="#8b949e",
                    alpha=0.18, lw=0, label=f"{base[0]}-{base[-1]} range")
    ax.plot(grid, np.median(bs, axis=0), color="white", ls="--", lw=1.8,
            label="baseline median")
    ax.set_xlabel("day of year", fontsize=10, color="white")
    ax.set_ylabel("cumulative share of that year's cuts", fontsize=10,
                  color="white")
    ax.set_ylim(0, 1)
    ax.set_title("Cumulative timing", fontsize=11.5, fontweight="bold",
                 color="white", loc="left")
    ax.legend(fontsize=8, frameon=True, facecolor="#161b22", edgecolor="none",
              labelcolor="white", loc="upper left")

    # A distribution-level summary that needs no threshold: the median cut date,
    # and how far the current year's whole distribution sits from the baseline.
    med_b = float(np.median([np.median(cuts.loc[cuts.year == y, "doy"])
                             for y in base]))
    med_c = float(np.median(cuts.loc[cuts.year == cur, "doy"]))
    q75_b = float(np.median([np.percentile(cuts.loc[cuts.year == y, "doy"], 75)
                             for y in base]))
    q75_c = float(np.percentile(cuts.loc[cuts.year == cur, "doy"], 75))
    ax.text(0.985, 0.06,
            f"median cut  DOY {med_b:.0f} \u2192 {med_c:.0f}  "
            f"({med_c - med_b:+.0f} d)\n"
            f"75th pct.   DOY {q75_b:.0f} \u2192 {q75_c:.0f}  "
            f"({q75_c - q75_b:+.0f} d)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=9,
            color="#c9d1d9",
            bbox=dict(facecolor="#161b22", edgecolor="#30363d",
                      boxstyle="round,pad=0.4"))

    _style(axes)
    n_parcels = d.parcel_id.nunique()
    fig.suptitle(f"{REGIME_TITLE[regime]}   ·   {site}   "
                 f"(n = {n_parcels} parcels)", fontsize=11.5,
                 fontweight="bold", color="#ff3b30")
    fig.supxlabel("Distributions, not counts: both panels are normalised within "
                  "the year, so they survive the fact that detection rate "
                  "differs between years.", fontsize=7, color="#6e7681")
    _save(fig, site, f"cut_timing_{regime}")
    plt.show()
    return cuts


def figure_cut_count(df, site, regime="intensive", cfg=None, n_rep=10, seed=0):
    """Distribution of the NUMBER of cuts per parcel, under matched observing
    conditions.

    READ THE CAVEAT BEFORE THE FIGURE. Cut count is the confounded metric:
    detection sensitivity rises from ~1.3 cuts at <15 clear images to ~2.9 at
    35-45, and drought summers are sunny, so the drought year has roughly twice
    the images of a cloudy baseline year. A raw count distribution therefore
    compares observing conditions as much as mowing.

    This figure avoids that by thinning every parcel-year to the SAME number of
    observations before detecting, repeated `n_rep` times. Both are shown: the
    raw distribution in outline, the matched one filled. Where they differ is
    the size of the artefact.

    Counts are small integers, so this is a stacked share rather than a violin —
    a kernel density over the values 0..5 would invent a smoothness that is not
    there.
    """
    C = cfg or CFG
    cur = C.project.current_year
    codes = C.landuse.intensive if regime == "intensive" else C.landuse.extensive
    d = df[df.lnf_code.isin(codes)] if "lnf_code" in df else df
    from .detect import season_end
    if C.detect.truncate_to_current:
        d = d[d.doy <= season_end(d, C)]
    d = d[d.doy >= C.detect.season[0]]

    raw = detect_all(d, C, verbose=False)

    rng = np.random.default_rng(seed)
    per = d.groupby(["parcel_id", "year"]).size()
    target = int(per.groupby(level=1).median().min())
    print(f"  thinning every parcel-year to {target} observations, "
          f"{n_rep} repeats")
    reps = []
    for _ in range(n_rep):
        keep = [g if len(g) <= target
                else g.sample(target, random_state=int(rng.integers(1 << 31)))
                for _, g in d.groupby(["parcel_id", "year"])]
        reps.append(detect_all(pd.concat(keep), C, verbose=False)
                    [["parcel_id", "year", "n_cuts"]])
    matched = pd.concat(reps)

    yrs = sorted(raw.year.unique())
    maxc = int(min(6, max(raw.n_cuts.max(), matched.n_cuts.max())))
    levels = list(range(0, maxc + 1))

    def shares(t):
        out = {}
        for y in yrs:
            v = t.loc[t.year == y, "n_cuts"].dropna()
            out[y] = [float((v == k).mean()) for k in levels[:-1]] + \
                     [float((v >= levels[-1]).mean())]
        return pd.DataFrame(out, index=[str(k) for k in levels[:-1]] +
                            [f"{levels[-1]}+"])

    sr, sm = shares(raw), shares(matched)

    plt.style.use("dark_background")
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2), dpi=180,
                             constrained_layout=True)
    ramp = plt.cm.viridis(np.linspace(0.15, 0.9, len(levels)))
    x = np.arange(len(yrs))

    for ax, tab, title in ((axes[0], sr, "Raw — confounded by image count"),
                           (axes[1], sm, f"Matched to {target} images per "
                                         "parcel-year")):
        bottom = np.zeros(len(yrs))
        for i, lab in enumerate(tab.index):
            vals = tab.loc[lab].to_numpy()
            ax.bar(x, vals, bottom=bottom, width=0.72, color=ramp[i],
                   edgecolor="none", label=lab)
            bottom += vals
        ax.set_xticks(x)
        ax.set_xticklabels(yrs)
        # Mark the current year on its LABEL rather than with a line through
        # the bar — a vertical rule over a stacked bar reads as a divider.
        for t, y in zip(ax.get_xticklabels(), yrs):
            if y == cur:
                t.set_color("#ff3b30"); t.set_fontweight("bold")
        ax.bar(x[[i for i, y in enumerate(yrs) if y == cur]], [0.018],
               bottom=[-0.028], width=0.72, color="#ff3b30",
               edgecolor="none", clip_on=False, zorder=5)
        ax.set_ylim(-0.03, 1)
        ax.set_ylabel("share of parcels", fontsize=10, color="white")
        ax.set_title(title, fontsize=11, fontweight="bold", color="white",
                     loc="left")

    axes[1].legend(title="cuts", fontsize=8, title_fontsize=8, frameon=True,
                   facecolor="#161b22", edgecolor="none", labelcolor="white",
                   loc="center left", bbox_to_anchor=(1.01, 0.5))
    mr = raw.groupby("year").n_cuts.mean()
    mm = matched.groupby("year").n_cuts.mean()
    axes[1].text(0.02, -0.16,
                 f"mean cuts, {cur}: raw {mr.get(cur, float('nan')):.2f}  ->  "
                 f"matched {mm.get(cur, float('nan')):.2f}   |   baseline "
                 f"{mr.drop(cur, errors='ignore').mean():.2f}  ->  "
                 f"{mm.drop(cur, errors='ignore').mean():.2f}",
                 transform=axes[1].transAxes, fontsize=8.5, color="#c9d1d9",
                 va="top")

    _style(axes)
    fig.suptitle(f"{REGIME_TITLE[regime]}   ·   {site}   ·   cuts per parcel "
                 f"({cur} marked in red)", fontsize=11.5, fontweight="bold",
                 color="#ff3b30")
    fig.supxlabel("Counts are small integers, so this is a stacked share, not a "
                  "violin: a kernel density over 0..5 would invent smoothness "
                  "that is not there.", fontsize=7, color="#6e7681")
    _save(fig, site, f"cut_count_{regime}")
    plt.show()
    return {"raw": sr, "matched": sm}


def figure_purity_sensitivity(site, regime="intensive", path=None, cfg=None,
                              min_px=200):
    """Does the anomaly depend on where the purity threshold is set?

    Reads the grouped export from `drought ee sensitivity`: mean index and
    pixel count per purity class per year. The mean at threshold t is the
    pixel-weighted combination of all classes at or above t — arithmetically
    identical to re-extracting at threshold t, but from one pass instead of one
    per threshold.

    The parcels that enter are those with at least one cell at gee.strict_purity.
    At the default 0.0 that is every parcel: the sweep covers the population the
    anomaly is reported for, and a drift may be contamination or composition. At
    0.90 the set is fixed large fields, and any drift is contamination alone.

    What to read: if the anomaly is flat from the strict end down to the chosen
    threshold, the choice is not a lever on the result — which is a stronger
    statement than any claim that 0.9 is optimal. Where it starts to bend is
    where contamination begins to matter.
    """
    C = cfg or CFG
    cur = C.project.current_year
    thr_used = C.mask.purity_min
    from .config import ROOT
    if path is None:
        cands = sorted((ROOT / C.paths.gee).glob(
            f"purity_sensitivity_{site}_{regime}*.csv"))
        if not cands:
            raise FileNotFoundError(
                f"no purity_sensitivity_{site}_{regime}*.csv in "
                f"{ROOT / C.paths.gee} — run `drought ee sensitivity` first")
        path = cands[0]
    from .provenance import verify
    verify(path, "sensitivity")            # refuses a stale export, naming why
    d = pd.read_csv(path)
    for c in ("year", "pclass", "purity_lo", "mean", "n_px"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d.dropna(subset=["mean", "n_px"])

    thresholds = sorted(d.purity_lo.unique())
    rows = []
    for t in thresholds:
        sub = d[d.purity_lo >= t - 1e-9]
        per_year = (sub.assign(w=sub["mean"] * sub.n_px)
                       .groupby("year")[["w", "n_px"]].sum())
        per_year["v"] = per_year.w / per_year.n_px
        if cur not in per_year.index:
            continue
        base = per_year.drop(index=cur)["v"]
        if base.empty:
            continue
        rows.append({"threshold": t, "current": per_year.loc[cur, "v"],
                     "baseline": base.median(),
                     "anomaly": per_year.loc[cur, "v"] - base.median(),
                     "n_px": per_year.loc[cur, "n_px"]})
    t = pd.DataFrame(rows)
    t = t[t.n_px >= min_px]
    if t.empty:
        print("not enough pixels at any threshold")
        return None

    ref = t.iloc[(t.threshold - thr_used).abs().argmin()]
    strict = t.iloc[-1]

    plt.style.use("dark_background")
    fig, axes = plt.subplots(2, 1, figsize=(7.6, 7.4), dpi=180, sharex=True,
                             constrained_layout=True,
                             gridspec_kw={"height_ratios": [2.2, 1]})

    ax = axes[0]
    ax.plot(t.threshold, t.anomaly, "o-", ms=4, lw=2.2, color="#ff3b30")
    ax.axvline(thr_used, color="white", ls="--", lw=1.6,
               label=f"threshold in use ({thr_used:.2f})")
    ax.axhline(strict.anomaly, color="#8b949e", ls=":", lw=1.4,
               label=f"strictest available ({strict.threshold:.2f})")
    ax.fill_between(t.threshold, strict.anomaly * 0.95, strict.anomaly * 1.05,
                    color="#8b949e", alpha=0.18, lw=0,
                    label="within 5% of strictest")
    ax.set_ylabel(f"{cur} anomaly, {C.detect.index}", fontsize=10, color="white")
    ax.set_title("Does the answer depend on the threshold?", fontsize=11.5,
                 fontweight="bold", color="white", loc="left")
    ax.legend(fontsize=8, frameon=True, facecolor="#161b22", edgecolor="none",
              labelcolor="white", loc="lower right")
    ax.text(0.03, 0.97,
            f"at {ref.threshold:.2f}:  {ref.anomaly:+.4f}\n"
            f"at {strict.threshold:.2f}:  {strict.anomaly:+.4f}\n"
            f"difference:  {100 * abs(ref.anomaly - strict.anomaly) / abs(strict.anomaly):.1f}% "
            "of the anomaly",
            transform=ax.transAxes, ha="left", va="top", fontsize=9,
            color="#c9d1d9",
            bbox=dict(facecolor="#161b22", edgecolor="#30363d",
                      boxstyle="round,pad=0.4"))

    ax = axes[1]
    ax.plot(t.threshold, t.n_px, "o-", ms=4, lw=2.0, color="#58a6ff")
    ax.axvline(thr_used, color="white", ls="--", lw=1.6)
    ax.set_yscale("log")
    ax.set_xlabel("purity threshold", fontsize=10, color="white")
    ax.set_ylabel("cells contributing", fontsize=10, color="#58a6ff")
    ax.tick_params(axis="y", colors="#58a6ff")
    ax.set_title("How much data is left, and how close to the true area",
                 fontsize=11, fontweight="bold", color="white", loc="left")

    # The area comparison lives in mask.area_table() / mask_figures.area_check() now: it needs the master
    # parcel set, which this module does not have.

    _style(axes)
    fig.suptitle(f"Purity threshold sensitivity  ·  {site}  ·  {regime}",
                 fontsize=11.5, fontweight="bold", color="#ff3b30")
    strict_p = C.gee.strict_purity
    who = ("Every parcel enters (gee.strict_purity 0): at high thresholds small fields "
           "drop out, so a drift can be contamination or composition."
           if not strict_p else
           f"Fixed parcel set (a cell >= {strict_p:g} required): the threshold changes "
           "which CELLS are used, never which FIELDS.")
    fig.supxlabel(f"{who}\nA flat curve means the choice is not a lever on the result.",
                  fontsize=7, color="#6e7681")
    _save(fig, site, f"purity_sensitivity_{regime}")
    plt.close(fig)

    print(f"{'threshold':>10}{'anomaly':>11}{'vs strictest':>14}{'cells':>10}")
    for _, r in t.iterrows():
        dv = 100 * abs(r.anomaly - strict.anomaly) / abs(strict.anomaly)
        mark = "  <- in use" if abs(r.threshold - thr_used) < 1e-9 else ""
        print(f"{r.threshold:>10.2f}{r.anomaly:>11.4f}{dv:>13.1f}%"
              f"{int(r.n_px):>10,}{mark}")
    return t
