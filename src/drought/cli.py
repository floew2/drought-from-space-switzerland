"""The `drought` command: one entry point for the whole analysis.

    drought demo [--fresh]                    everything, on synthetic data
    drought run-all [--synthetic]             every step that has input data
    drought synthetic generate [--seed N] [--scale X]
    drought mask build | check | figures
    drought trajectory plot [--site S] [--index NDMI]
    drought precip [--download]
    drought cuts detect --site S
    drought cuts score --site S
    drought gallery
    drought notebook
    drought ee check|stamp-mask|trajectory|sensitivity|cuts|run-all
    drought sensitivity plot [--site S] [--regime R]
    drought provenance register|verify|show FILES [--kind KIND]
    drought results show
    drought config show [KEY]

Where settings come from, in increasing precedence:

1. `config.yaml`: every default, in one place.
2. `--config FILE`: a whole alternative configuration.
3. `--synthetic`: switches every data and output path to the synthetic profile.
4. `--set KEY=VALUE` (repeatable): one setting, e.g.
   `--set mask.purity_min=0.8`. Only existing keys are accepted.

The notebook does not go through this: it passes explicit arguments to the
library functions. The libraries themselves (config, io, clean, metrics, …)
have no command line.

Overrides are applied before any analysis module is imported. Some modules
copy settings into constants at import, and the configuration refuses a change
that could no longer reach them, rather than letting it silently do nothing.
"""

from __future__ import annotations

import argparse
import sys

import pandas as pd
import yaml


# --------------------------------------------------------------------------- #
# Commands. Each imports the analysis modules it needs only when it runs,
# so overrides are in place first.
# --------------------------------------------------------------------------- #
def _demo(a) -> None:
    from . import demo
    demo.run(fresh=a.fresh)


def _synthetic_generate(a) -> None:
    from . import synthetic
    synthetic.generate(seed=a.seed, scale=a.scale)


def _mask_build(a) -> None:
    from . import config, mask, results
    master = mask.build_mask(shp_dir=a.shp_dir)
    results.update("mask", {"funnel": master.attrs["funnel"], "response": mask.psf_params()},
                   inputs=[config.mask_path()])


def _mask_check(a) -> None:
    from . import config, mask, results
    master = mask.master_parcels(a.shp_dir, verbose=False)
    table = mask.area_table(master, site=a.site)
    table.pop("_window")
    check = mask.psf_model_check(master)
    results.update("mask_area", {
        "window": a.site or "full extent",
        "categories": {n: {"parcels": d["parcels"], "polygon_ha": round(d["polygon_ha"], 1),
                           "ratio": {f"{t:g}": round(r, 4) for t, r in d["ratio"].items()}}
                       for n, d in table.items()},
        "psf_model_check": {f"{c}@{t:g}": {k: round(float(v), 4) for k, v in row.items()}
                            for (c, t), row in check.iterrows()}},
        inputs=[config.mask_path()])


def _mask_overview(a) -> None:
    from . import mask, mask_figures
    master = mask.master_parcels(a.shp_dir, verbose=False)
    mask_figures.figure_mask_overview(thresholds=tuple(a.thresholds), master=master)


def _mask_figures(a) -> None:
    from . import mask, mask_figures
    master = mask.master_parcels(a.shp_dir, verbose=False)
    areas = {n: float(master[master["_code"].isin(c)].geometry.area.sum() / 1e4)
             for n, c in mask.GROUPS.items()}
    mask_figures.figure_psf()
    mask_figures.figure_parcel_size(master=master)
    mask_figures.area_check(master=master)
    for site in a.sites:
        mask_figures.figure_purity_explainer(site=site, areas=areas)
        mask_figures.figure_purity_map(site, group="intensive")


def _trajectory_plot(a) -> None:
    from . import results, trajectory
    src = a.path or trajectory._find_csv(a.site)
    res = trajectory.plot(site=a.site, path=src, index=a.index)
    if res is not None:
        print(res[["site", "regime", "integral_pct", "rank", "upto", "n_days"]]
              .to_string(index=False))
        results.update("condition", {
            f"{r.site}/{r.regime}": {"index": trajectory.INDEX,
                                     "integral_pct": round(float(r.integral_pct), 1),
                                     "rank": int(r.rank), "n_baseline_years": int(r.n_baseline_years),
                                     "upto_doy": int(r.upto), "shared_days": int(r.n_days)}
            for r in res.itertuples()}, inputs=[src], merge=True)


def _precip(a) -> None:
    from . import config, precip
    if a.download:
        precip.main()
        return
    try:
        daily, abbr, lat = precip.load_local(config.station_dir())
    except FileNotFoundError:
        sys.exit(f"no station files in {config.station_dir()}; "
                 "use --download to fetch them from MeteoSwiss")
    from . import results
    label = "Synthetic station" if config.active_profile() == "synthetic" else None
    res = precip.main(daily=daily, abbr=abbr, lat=lat, station_label=label)
    res.pop("file", None)
    results.update("precipitation", res,
                   inputs=sorted(config.station_dir().glob("ogd-smn_*_d_*.csv")))


def _cuts_detect(a) -> None:
    import drought as dr
    df = dr.load_site(a.site, labelled=a.labelled)
    ev = dr.detect.detect_all(df)
    print(dr.metrics.summarise(ev).to_string())


def _cuts_score(a) -> None:
    import drought as dr

    from . import changepoint, evaluate, labels
    df = dr.load_site(a.site, verbose=False, labelled=a.labelled)
    lab = labels.check_labels(labels.load_labels(a.site), df)
    evaluate.assert_wired()
    print("heuristic (frozen configuration):")
    evaluate.score(lab, df)
    print(evaluate.score_by_year(lab, df)[["n_parcel_years", "precision", "recall"]]
          .to_string())
    cp = changepoint.score_detector(lab, df, changepoint.detect_changepoints)
    print(f"\nchangepoint alternative: precision {cp['precision']}  recall "
          f"{cp['recall']}  F1 {cp['f1']}")
    from . import config, results
    h = evaluate.score(lab, df, verbose=False)
    by_year = evaluate.score_by_year(lab, df)
    results.update(f"detector.{a.site}", {
        "labelled_parcel_years": len(lab), "precision": h["precision"],
        "recall": h["recall"], "f1": h["f1"],
        "recall_per_year": [float(by_year.recall.min()), float(by_year.recall.max())],
        "changepoint_f1": cp["f1"]},
        inputs=sorted((config.gee_dir(a.site) / ("labelled" if a.labelled else ""))
                      .glob("cuts_*.csv")))


def _cuts_errors(a) -> None:
    import drought as dr

    from . import evaluate, labels
    df = dr.load_site(a.site, verbose=False, labelled=a.labelled)
    lab = labels.check_labels(labels.load_labels(a.site), df)
    by_year = evaluate.score_by_year(lab, df)
    cur = int(dr.CFG.project.current_year)
    rest = by_year.drop(index=cur, errors="ignore")[["tp", "fp", "fn"]].sum()
    p = rest.tp / (rest.tp + rest.fp)
    r = rest.tp / (rest.tp + rest.fn)
    print(f"\nwithout {cur}: precision {p:.3f}  recall {r:.3f}  F1 {2 * p * r / (p + r):.3f}")
    if cur in by_year.index:
        y = by_year.loc[cur]
        print(f"{cur} alone: {int(y.tp)} hits, {int(y.fp)} false positives, {int(y.fn)} misses")
    fp = evaluate.false_positive_positions(lab, df)
    bins = list(range(60, 320, 20))
    if not fp.empty:
        fp["window"] = pd.cut(fp.doy, bins=bins, right=False)
        tab = fp.pivot_table(index="window", columns=fp.year == cur, values="doy",
                             aggfunc="count", fill_value=0, observed=False)
        tab.columns = [str(cur) if c else "other years" for c in tab.columns]
        print("\nfalse positives by day of year:")
        print(tab.to_string())
    miss = evaluate.miss_positions(lab, df)
    if not miss.empty:
        print(f"\nmisses: {len(miss)}, of them last cuts of their season: "
              f"{int(miss.is_last.sum())}")


def _cuts_labels(a) -> None:
    from . import labels, metrics
    lab = labels.load_labels(a.site)
    t = metrics.labels_last_cut_paired(lab, end=a.end)
    if t.empty:
        sys.exit("no parcel labelled in the current year and a baseline year")
    print(f"hand labels only, every year cut at DOY {a.end}:")
    print(f"  {len(t)} parcels labelled in 2026 and >= 1 baseline year "
          f"(median {t.n_baseline_years.median():.0f} baseline years each)")
    print(f"  baseline mean last cut  DOY {t.baseline.mean():.1f}")
    print(f"  2026 mean last cut      DOY {t.current.mean():.1f}")
    print(f"  median shift            {t.shift_days.median():+.1f} days")
    print(f"  parcels stopping sooner {int((t.shift_days < 0).sum())} of {len(t)}")


def _notebook(a) -> None:
    from . import notebook
    notebook.run(check=a.check)


def _gallery(a) -> None:
    from . import gallery
    gallery.build()


def _run_all(a) -> None:
    from . import config
    if a.synthetic:
        _demo(argparse.Namespace(fresh=False))
        return
    steps = [
        ("mask", lambda: _mask_build(argparse.Namespace(shp_dir=None))),
        ("trajectory", lambda: _trajectory_plot(
            argparse.Namespace(site=None, path=None, index=None))),
        ("precipitation", lambda: _precip(argparse.Namespace(download=False))),
    ]
    for site in vars(config.CFG.sites):
        if config.labels_path(site).exists():
            steps.append((f"cuts {site}",
                          lambda s=site: _cuts_score(argparse.Namespace(site=s))))
    failed = []
    for name, step in steps:
        print(f"\n=== {name} " + "=" * (60 - len(name)))
        try:
            step()
        except (FileNotFoundError, SystemExit) as exc:
            print(f"  [skipped] {exc}")
            failed.append(name)
    if failed:
        print(f"\nskipped for missing input: {', '.join(failed)}")


def _sensitivity_plot(a) -> None:
    from . import figures, results
    from .config import CFG

    sites = a.site or list(CFG.gee.sensitivity.windows.to_dict())
    regimes = a.regime or list(CFG.analysis.regimes)
    out, used = {}, CFG.mask.purity_min
    for site in sites:
        for regime in regimes:
            print(f"\n{site} / {regime}")
            try:
                t = figures.figure_purity_sensitivity(site, regime)
            except FileNotFoundError as exc:
                print(f"  [skipped] {exc}")
                continue
            if t is None:
                continue
            strict = t.iloc[-1]
            ref = t.iloc[(t.threshold - used).abs().argmin()]
            hi = t[t.threshold >= 0.5]
            spread = float((hi.anomaly.max() - hi.anomaly.min()) / abs(strict.anomaly))
            out[f"{site}/{regime}"] = {
                "anomaly_at_threshold": round(float(ref.anomaly), 4),
                "threshold": float(ref.threshold),
                "anomaly_strictest": round(float(strict.anomaly), 4),
                "strictest": float(strict.threshold),
                "spread_0p5_up_pct_of_anomaly": round(100 * spread, 1)}
    if out:
        results.update("sensitivity", out)
        print("\n  spread of the anomaly over thresholds 0.5-0.95, as % of the anomaly:")
        for k, v in out.items():
            print(f"    {k:<24}{v['spread_0p5_up_pct_of_anomaly']:>6.1f} %")


def _ee(a) -> None:
    from .config import CFG
    from .engine import cuts, imagery, sensitivity, session, stamp, tasks, trajectory

    if a.action == "check":
        print(f"  project {session.project()}; exports go to Drive folder "
              f"{CFG.gee.drive_folder!r}")
        session.init()
        session.drive_find("drought-from-space connectivity check")
        print("  Earth Engine and Google Drive are reachable")
        sha = imagery.check_mask_asset()
        print(f"  mask asset is stamped as the local mask ({sha}); CONFIG_HASH "
              f"{tasks.config_hash()}")
        return
    if a.action == "stamp-mask":
        stamp.stamp_mask()
        return
    if a.action == "rgb":
        from .engine import rgb
        rgb.run(sites=a.site)
        return
    imagery.check_mask_asset()                       # before building anything
    if a.action == "trajectory":
        exps = trajectory.exports(smoke=a.smoke)
    elif a.action == "sensitivity":
        exps = sensitivity.exports(sites=a.site, regimes=a.regime)
    elif a.action == "cuts":
        exps = cuts.exports(sites=a.site, only_years=a.year, smoke=a.smoke,
                            labelled_only=a.labelled_only)
    else:                                            # run-all
        exps = (trajectory.exports() + sensitivity.exports() + cuts.exports())
    paths = tasks.run(exps, dry_run=a.dry_run)
    if paths:
        n_reg = sum(1 for e in exps if e.kind and e.local in paths)
        print(f"  {len(paths)} file(s) downloaded, {n_reg} registered")


def _provenance(a) -> None:
    import json

    from . import provenance
    try:
        if a.action == "register":
            for f in a.files:
                print(f"  registered {f} -> {provenance.register(f, a.kind).name}")
        elif a.action == "verify":
            for f in a.files:
                provenance.verify(f, a.kind)
                print(f"  ok  {f}")
        else:
            for f in a.files:
                print(json.dumps(json.loads(provenance.sidecar_path(f).read_text()), indent=1))
    except provenance.ProvenanceError as exc:
        sys.exit(f"drought: {exc}")


def _results_show(a) -> None:
    import json

    from . import results
    try:
        doc = results.load(strict=not a.allow_stale)
    except FileNotFoundError:
        sys.exit(f"no results yet: {results.path()}")
    except RuntimeError as exc:
        sys.exit(f"drought: {exc}")
    print(f"# {results.path()}")
    print(json.dumps({k: v["data"] for k, v in doc["sections"].items()}, indent=1))


def _config_show(a) -> None:
    from . import config
    cfg = config.CFG.to_dict()
    if a.key:
        node = config.CFG.get(a.key, default=KeyError)
        if node is KeyError:
            sys.exit(f"no setting {a.key!r}")
        cfg = node.to_dict() if isinstance(node, config.Config) else node
    print(f"# profile: {config.active_profile()}")
    if isinstance(cfg, dict):
        print(yaml.safe_dump(cfg, sort_keys=False, allow_unicode=True).rstrip())
    else:
        print(cfg)


# --------------------------------------------------------------------------- #
# Parser
# --------------------------------------------------------------------------- #
def _parse_set(item: str) -> tuple[str, object]:
    if "=" not in item:
        raise argparse.ArgumentTypeError(f"expected KEY=VALUE, got {item!r}")
    key, raw = item.split("=", 1)
    return key.strip(), yaml.safe_load(raw)        # 0.8 -> float, [1, 2] -> list


def build_parser() -> argparse.ArgumentParser:
    """The argument parser for the `drought` command."""
    p = argparse.ArgumentParser(
        prog="drought", description="Reading a drought from space: the analysis pipeline.",
        epilog="Precedence: config.yaml < --config < --synthetic < --set.")
    p.add_argument("--config", metavar="FILE", help="use this YAML instead of config.yaml")
    p.add_argument("--synthetic", action="store_true",
                   help="use the synthetic data profile (data/synthetic, output/synthetic)")
    p.add_argument("--set", dest="overrides", metavar="KEY=VALUE", action="append",
                   type=_parse_set, default=[],
                   help="override one setting, e.g. --set mask.purity_min=0.8 (repeatable)")
    sub = p.add_subparsers(dest="command", required=True, metavar="COMMAND")

    c = sub.add_parser("demo", help="the whole pipeline on synthetic data, self-checking")
    c.add_argument("--fresh", action="store_true", help="regenerate the synthetic inputs")
    c.set_defaults(func=_demo)

    c = sub.add_parser("run-all", help="every step that has input data")
    c.add_argument("--synthetic", action="store_true", dest="synthetic_run",
                   help="same as `drought demo`")
    c.set_defaults(func=lambda a: _run_all(
        argparse.Namespace(synthetic=a.synthetic or a.synthetic_run)))

    syn = sub.add_parser("synthetic", help="synthetic inputs").add_subparsers(
        dest="action", required=True, metavar="ACTION")
    c = syn.add_parser("generate", help="write every synthetic input")
    c.add_argument("--seed", type=int, help="master seed (default: config)")
    c.add_argument("--scale", type=float, help="share of the real parcel count (default: config)")
    c.set_defaults(func=_synthetic_generate)

    m = sub.add_parser("mask", help="the grassland purity mask").add_subparsers(
        dest="action", required=True, metavar="ACTION")
    for name, func, text in (("build", _mask_build, "build the four-band purity raster"),
                             ("check", _mask_check, "area cost per threshold; PSF model check"),
                             ("figures", _mask_figures, "the mask figures"),
                             ("overview", _mask_overview,
                              "imagery, field outlines and purity, one row per site")):
        c = m.add_parser(name, help=text)
        c.add_argument("--shp-dir", help="folder with the yearly shapefiles (default: config)")
        if name == "check":
            c.add_argument("--site", help="restrict to one site window")
        if name == "overview":
            c.add_argument("--thresholds", nargs="+", type=float, default=[0.0, 0.5, 0.9],
                           help="purity thresholds, one column each (default: 0 0.5 0.9)")
        if name == "figures":
            c.add_argument("--sites", nargs="+", default=["linthebene"],
                           help="sites for the per-site figures")
        c.set_defaults(func=func)

    t = sub.add_parser("trajectory", help="site-mean index trajectories").add_subparsers(
        dest="action", required=True, metavar="ACTION")
    c = t.add_parser("plot", help="trajectory and season-integral figures")
    c.add_argument("--site", help="one site (default: all in the export)")
    c.add_argument("--index", help="NDMI, NDVI, kNDVI or NDRE (default: config detect.index)")
    c.add_argument("--path", help="trajectory CSV (default: found under paths.gee)")
    c.set_defaults(func=_trajectory_plot)

    se = sub.add_parser("sensitivity", help="purity-threshold sensitivity").add_subparsers(
        dest="action", required=True, metavar="ACTION")
    c = se.add_parser("plot", help="anomaly at every purity threshold, per site and system")
    c.add_argument("--site", action="append", help="repeatable; default: all")
    c.add_argument("--regime", action="append", help="repeatable; default: all")
    c.set_defaults(func=_sensitivity_plot)

    c = sub.add_parser("precip", help="cumulative precipitation and P - ET0")
    c.add_argument("--download", action="store_true",
                   help="fetch from MeteoSwiss instead of reading paths.station")
    c.set_defaults(func=_precip)

    cu = sub.add_parser("cuts", help="cut detection").add_subparsers(
        dest="action", required=True, metavar="ACTION")
    for name, func, text in (("detect", _cuts_detect, "detect cuts, per-year summary"),
                             ("score", _cuts_score, "score against the hand labels"),
                             ("errors", _cuts_errors,
                              "where the false positives and misses fall")):
        c = cu.add_parser(name, help=text)
        c.add_argument("--site", required=True)
        c.add_argument("--labelled", action="store_true",
                       help="read data/gee/<site>/labelled/ (from `ee cuts --labelled-only`)")
        c.set_defaults(func=func)

    c = cu.add_parser("labels", help="last-cut shift from the hand labels alone (no detector)")
    c.add_argument("--site", required=True)
    c.add_argument("--end", type=int, default=228, help="comparison end DOY (default 228)")
    c.set_defaults(func=_cuts_labels)

    sub.add_parser("gallery", help="refresh docs/figures/model").set_defaults(func=_gallery)
    nb = sub.add_parser("notebook", help="run notebooks/walkthrough.ipynb on synthetic data")
    nb.add_argument("--check", action="store_true", help="run a copy; fail on any error (CI)")
    nb.set_defaults(func=_notebook)


    e = sub.add_parser("ee", help="run the Earth Engine exports from Python").add_subparsers(
        dest="action", required=True, metavar="ACTION")
    e.add_parser("check", help="project, Drive and mask asset in order?").set_defaults(func=_ee)
    e.add_parser("stamp-mask", help="compare the uploaded asset with the local mask; "
                 "stamp it with the mask's hash if they agree").set_defaults(func=_ee)
    c = e.add_parser("rgb", help="false-colour composite per site -> data/gee/<site>/rgb_<site>.tif "
                     "(for `drought mask overview`)")
    c.add_argument("--site", action="append", help="repeatable; default: all")
    c.set_defaults(func=_ee)
    for name, text in (("trajectory", "site-mean index series -> data/gee/"),
                       ("sensitivity", "purity-class series per site and system -> data/gee/"),
                       ("cuts", "per-parcel series per year -> data/gee/<site>/"),
                       ("run-all", "trajectory, sensitivity and cuts, one after another")):
        c = e.add_parser(name, help=text)
        if name in ("sensitivity", "cuts"):
            c.add_argument("--site", action="append", help="repeatable; default: all")
        if name == "sensitivity":
            c.add_argument("--regime", action="append", help="repeatable; default: all")
        if name == "cuts":
            c.add_argument("--year", action="append", type=int, help="repeatable")
        if name in ("trajectory", "cuts"):
            c.add_argument("--smoke", action="store_true",
                           help="a small run to measure cost; not registered")
        if name == "cuts":
            c.add_argument("--labelled-only", action="store_true",
                           help="only the hand-labelled parcel-years: enough to re-score")
        c.add_argument("--dry-run", action="store_true",
                       help="check the mask asset and list the exports, start nothing")
        c.set_defaults(func=_ee)

    pv = sub.add_parser("provenance", help="provenance records of data files").add_subparsers(
        dest="action", required=True, metavar="ACTION")
    for name, text in (("register", "record an Earth Engine export on arrival, after checking it"),
                       ("verify", "check files against the current configuration and mask"),
                       ("show", "print a file's record")):
        c = pv.add_parser(name, help=text)
        c.add_argument("files", nargs="+")
        if name != "show":
            c.add_argument("--kind", required=True, choices=["trajectory", "cuts", "sensitivity", "station", "mask"])
        c.set_defaults(func=_provenance)

    r = sub.add_parser("results", help="output/results.json").add_subparsers(
        dest="action", required=True, metavar="ACTION")
    c = r.add_parser("show", help="print every headline number, refusing stale sections")
    c.add_argument("--allow-stale", action="store_true")
    c.set_defaults(func=_results_show)

    c = sub.add_parser("config", help="inspect the effective configuration").add_subparsers(
        dest="action", required=True, metavar="ACTION").add_parser(
        "show", help="print it, after --config, --synthetic and --set")
    c.add_argument("key", nargs="?", help="dotted key, e.g. mask.purity_min")
    c.set_defaults(func=_config_show)
    return p


def main(argv: list[str] | None = None) -> None:
    """Entry point of the `drought` console script."""
    args = build_parser().parse_args(argv)
    import matplotlib
    matplotlib.use("Agg")                         # a command line never opens windows

    from . import config
    try:
        if args.config:
            config.use_config(args.config)
        synthetic = args.synthetic or args.command == "demo" or \
            getattr(args, "synthetic_run", False)
        if synthetic:
            config.activate("synthetic")
        for key, value in args.overrides:
            config.override(key, value)
    except KeyError as exc:
        sys.exit(f"drought: {exc.args[0]}")
    if synthetic:
        from .synthetic import stamp_figures
        stamp_figures()
    from .provenance import ProvenanceError
    try:
        args.func(args)
    except ProvenanceError as exc:
        sys.exit(f"drought: refused — {exc}")


if __name__ == "__main__":
    main()
