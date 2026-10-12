"""`make demo`: the whole pipeline, end to end, on synthetic data.

    python -m drought.demo            # generate (if needed) and run
    python -m drought.demo --fresh    # regenerate the synthetic inputs first

Runs the same code the real analysis uses — mask build, area accounting,
trajectories and season integrals, precipitation, cut detection and scoring —
on synthetic inputs calibrated to the real distributions. Outputs go to
output/synthetic/; every figure is stamped SYNTHETIC.

The closing table puts each synthetic number beside the real one it was
calibrated to (or, for emergent quantities, should resemble). Calibrated rows
reproducing their target shows the pipeline is wired correctly; it is not new
evidence about the drought.
"""

from __future__ import annotations

import argparse
import json
import time

import numpy as np
import pandas as pd


def _funnel(master, truth: pd.DataFrame) -> dict:
    """The master-set funnel, against what the generator actually planted."""
    kept = set(master["objectid"].astype(int))
    t = truth.assign(kept=truth.objectid.astype(int).isin(kept))
    by = t.groupby("status").kept.agg(["size", "sum"])
    stable = by.loc["stable"] if "stable" in by.index else pd.Series({"size": 0, "sum": 0})
    lost_stable = t[(t.status == "stable") & ~t.kept]
    return {
        "grassland_2025": int(len(t)),
        "remain": int(t.kept.sum()),
        "remain_share": round(float(t.kept.mean()), 3),
        "planted": {k: int(v) for k, v in t.status.value_counts().items()},
        "stable_dropped": int(stable["size"] - stable["sum"]),
        "stable_dropped_max_ha": round(float(lost_stable.area_ha.max()), 3)
        if len(lost_stable) else None,
    }


def _geometry(data: dict) -> dict:
    out = {}
    for n, (a, shp) in data.items():
        out[n] = {"parcels": int(len(a)), "median_ha": round(float(np.median(a)), 3),
                  "lt_05": round(float((a < 0.5).mean()), 3),
                  "median_shape": round(float(np.median(shp)), 2),
                  "shape_gt10": round(float((shp > 10).mean()), 4),
                  "shape_gt16": int((shp > 16).sum())}
    return out


def _cuts(site: str) -> dict:
    import drought as dr
    from drought import changepoint, evaluate, figures, labels

    df = dr.load_site(site, verbose=False)
    lab = labels.check_labels(labels.load_labels(site), df)
    h = evaluate.score(lab, df, verbose=False)
    by_year = evaluate.score_by_year(lab, df)
    cp = changepoint.score_detector(lab, df, changepoint.detect_changepoints)
    figures.figure_last_cut(df, site, regime="intensive")
    return {"labelled_parcel_years": len(lab), "precision": h["precision"],
            "recall": h["recall"], "f1": h["f1"],
            "recall_per_year": (float(by_year.recall.min()), float(by_year.recall.max())),
            "changepoint_f1": cp["f1"]}


def run(fresh: bool = False) -> dict:
    """Generate (if needed) and run everything; return and write the summary."""
    import matplotlib

    matplotlib.use("Agg")
    from drought import config, synthetic
    from drought.synthetic import targets as T

    t0 = time.time()
    config.activate("synthetic")
    synthetic.stamp_figures()
    if fresh or not synthetic.exists():
        gen = synthetic.generate()
    else:
        gen = {"note": "existing synthetic data reused; --fresh to regenerate"}
        print("synthetic data present — reusing (python -m drought.demo --fresh to regenerate)")

    from drought import mask, mask_figures, precip, trajectory

    summary = {"generator": gen}

    print("\n=== 1  grassland mask ==========================================")
    master = mask.build_mask()
    truth = pd.read_csv(config.parcels_dir() / "parcels_truth.csv")
    summary["funnel"] = _funnel(master, truth)
    area = mask_figures.area_check(master=master)
    summary["area_ratio"] = {n: (round(v["ratio"][0.0], 2), round(v["ratio"][0.9], 2))
                             for n, v in area.items()}
    summary["psf_model_check"] = {
        f"{c}@{t:g}": {k: round(float(v), 4) for k, v in row.items()}
        for (c, t), row in mask.psf_model_check(master).iterrows()}
    summary["psf"] = mask_figures.figure_psf()
    summary["geometry"] = _geometry(mask_figures.figure_parcel_size(master=master))
    areas = {n: float(master[master["_code"].isin(c)].geometry.area.sum() / 1e4)
             for n, c in mask.GROUPS.items()}
    mask_figures.figure_purity_explainer(site="linthebene", group="intensive", areas=areas)
    mask_figures.figure_purity_map("linthebene", group="intensive")

    print("\n=== 2  grassland condition =====================================")
    traj = trajectory.plot(path=config.ROOT / config.CFG.paths.gee
                           / "grassland_trajectory_multi.csv", index="NDMI")
    summary["trajectory"] = {
        f"{r.site}/{r.regime}": {"integral_pct": round(float(r.integral_pct), 1),
                                 "rank": int(r.rank), "upto": int(r.upto),
                                 "n_days": int(r.n_days)}
        for r in traj.itertuples()}

    print("\n=== 3  precipitation ===========================================")
    daily, abbr, lat = precip.load_local(config.station_dir())
    summary["precip"] = precip.main(daily=daily, abbr=abbr, lat=lat,
                                    station_label="Synthetic station")

    print("\n=== 4  cut detection ===========================================")
    summary["detector"] = {s: _cuts(s) for s in ("linthebene", "rheintal")}

    summary["runtime_s"] = round(time.time() - t0, 1)
    summary["failed_checks"] = check(summary, T)
    out = config.tables_dir() / "demo_summary.json"
    out.write_text(json.dumps(summary, indent=1, default=str))
    _report(summary, T)
    print(f"\n  summary -> {out.relative_to(config.ROOT)}")
    print(f"  figures -> {config.figures_dir().relative_to(config.ROOT)}/")
    print(f"  {summary['runtime_s']:.0f} s")
    if summary["failed_checks"]:
        raise SystemExit("demo checks FAILED:\n  " + "\n  ".join(summary["failed_checks"]))
    print("  all demo checks passed")
    return summary


def check(s: dict, T) -> list[str]:
    """What must hold for the demo to count as a pass. Returns the failures.

    Calibrated quantities must reproduce their targets; the physical invariant
    must hold; the response model must match the stated derivation. CI runs the
    demo, so a change that breaks any of these fails the build.
    """
    bad = []
    fn = s["funnel"]
    if fn["remain"] != fn["planted"].get("stable", 0) + fn["planted"].get("exempt_930", 0):
        bad.append(f"master set kept {fn['remain']} parcels, planted stable + exempt "
                   f"{fn['planted']}")
    if s["funnel"]["stable_dropped"]:
        bad.append(f"{s['funnel']['stable_dropped']} parcels planted as stable failed the "
                   "persistence test: the funnel must recover exactly what was planted")
    for (site, reg), v in T.INTEGRAL_PCT.items():
        t = s["trajectory"].get(f"{site}/{reg}")
        if t is None or abs(t["integral_pct"] - v) > 0.5 or t["rank"] != 1:
            bad.append(f"integral {site}/{reg}: {t} against {v}% rank 1")
    p = s["precip"]
    if abs(p["precip_mm"] - T.PRECIP["p_2026_mm"]) > 1.5 or p["precip_rank"] != 1:
        bad.append(f"precipitation {p['precip_mm']:.1f} mm rank {p['precip_rank']}")
    if abs(p["wb_mm"] - T.PRECIP["wb_2026_mm"]) > 1.5:
        bad.append(f"P - ET0 {p['wb_mm']:.1f} mm")
    for n in ("intensive", "extensive", "grazed"):
        g, r = s["geometry"][n], T.GEOMETRY[n]
        if abs(g["median_ha"] / r["median_ha"] - 1) > 0.10 or abs(g["lt_05"] - r["lt_05"]) > 0.04 \
                or abs(g["median_shape"] - r["median_shape"]) > 0.15:
            bad.append(f"geometry {n}: {g} against {r}")
    for n, (r0, _) in s["area_ratio"].items():
        if not r0 > 1.0:
            bad.append(f"area ratio at >0 for {n} is {r0}: the PSF must reach past every boundary")
    if round(s["psf"]["sigma"], 2) != 4.35 or round(s["psf"]["fwhm"], 1) != 12.7:
        bad.append(f"response model: {s['psf']}")
    return bad


def _report(s: dict, T) -> None:
    """Synthetic beside real, so a reader sees what the demo does and does not show."""
    rows = []

    def add(group, what, syn, real, kind):
        rows.append((group, what, syn, real, kind))

    g = s["geometry"]
    for n in ("intensive", "extensive", "mown", "grazed"):
        r = T.GEOMETRY[n]
        add("geometry", f"{n} median ha / <0.5 ha / shape",
            f"{g[n]['median_ha']:.2f} / {g[n]['lt_05']:.0%} / {g[n]['median_shape']:.2f}",
            f"{r['median_ha']:.2f} / {r['lt_05']:.0%} / {r['median_shape']:.2f}",
            "calibrated" if n in ("intensive", "extensive", "grazed") else "emergent")
    for n, (r0, r9) in T.AREA_RATIO.items():
        a0, a9 = s["area_ratio"][n]
        add("mask", f"{n} area ratio >0 / >=0.9", f"{a0:.2f} / {a9:.2f}",
            f"{r0:.2f} / {r9:.2f}", "emergent")
    f = s["funnel"]
    add("mask", "master set retained", f"{f['remain_share']:.1%}",
        f"{T.FUNNEL['remain'] / T.FUNNEL['grassland_2025']:.1%} (old filter)", "emergent")
    for (site, reg), v in T.INTEGRAL_PCT.items():
        t = s["trajectory"].get(f"{site}/{reg}")
        if t:
            add("condition", f"{site} {reg} integral (rank)",
                f"{t['integral_pct']:+.0f}% ({t['rank']}/8)", f"{v:+.0f}% (1/8)", "calibrated")
    p = s["precip"]
    add("precip", "P to 10 Aug, % of median (rank)",
        f"{p['precip_mm']:.0f} mm, {p['precip_pct_of_median']:.0f}% ({p['precip_rank']}/31)",
        "501 mm, 61% (1/31)", "calibrated")
    add("precip", "P - ET0 to 10 Aug", f"{p['wb_mm']:+.0f} mm", "-84 mm", "calibrated")
    for site, d in s["detector"].items():
        add("detector", f"{site} P / R / F1",
            f"{d['precision']:.2f} / {d['recall']:.2f} / {d['f1']:.2f}",
            f"{T.DETECTOR['precision']:.2f} / {T.DETECTOR['recall']:.2f} / "
            f"{T.DETECTOR['f1']:.2f}*", "emergent")
        add("detector", f"{site} changepoint F1", f"{d['changepoint_f1']:.2f}",
            f"{T.DETECTOR['changepoint_f1']:.2f}*", "emergent")

    print("\n=== synthetic against real ======================================")
    print(f"  {'':<10}{'quantity':<38}{'synthetic':<26}{'real':<22}")
    for grp, what, syn, real, kind in rows:
        print(f"  {grp:<10}{what:<38}{syn:<26}{real:<22}{kind}")
    print("\n  calibrated  the generator was solved to reproduce the real value")
    print("  emergent    not targeted; produced by the method acting on synthetic input")
    print("  *           real detector scores are on the previous 20 m mask, pooled sites")
    print(f"\n  master set: {f['remain']:,} of {f['grassland_2025']:,} kept; planted "
          + ", ".join(f"{k} {v}" for k, v in f["planted"].items())
          + ".")
    if f["stable_dropped"]:
        print(f"  {f['stable_dropped']} parcels planted as stable were dropped anyway "
              f"(largest {f['stable_dropped_max_ha']} ha).")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fresh", action="store_true", help="regenerate the synthetic inputs")
    run(fresh=ap.parse_args(argv).fresh)


if __name__ == "__main__":
    main()
