"""The `drought` command.

Pipeline:
    drought run [STAGE] [--aoi A]... [--force] [--dry-run]
                                    STAGE and everything it needs that is not
                                    done (default STAGE: report, i.e. everything)
Stages, one at a time (each fails with the command to run if an input is missing):
    drought fetch --aoi A [--years 2019-2026] [--no-vhi] [--dry-run]
    drought mask build | check | figures | overview
    drought condition | sensitivity | anomaly | cuts | vhi  --aoi A
    drought evaluate --aoi A
    drought precip [--download]
    drought report [--aoi A]... [--only PART]...
Looking closer:
    drought inspect detect | errors | field | labels --aoi A
Synthetic data and housekeeping:
    drought demo [--fresh]          the pipeline on synthetic data, self-checking
    drought synthetic generate      synthetic inputs only
    drought notebook [--check]      run notebooks/walkthrough.ipynb on synthetic data
    drought gallery                 refresh docs/figures/model
    drought provenance verify | show FILES
    drought results show            every headline number, refusing stale ones
    drought config show [KEY]

An AOI (`--aoi`) is a site from config.yaml or a GeoJSON file of any polygon
in Switzerland. Without `--aoi`, each stage runs for its configured sites.

Where settings come from, in increasing precedence: config.yaml; `--config
FILE`; `--synthetic` (every data and output path switches to the synthetic
profile); `--set KEY=VALUE` (repeatable, existing keys only). Overrides are
applied before any analysis module is imported.
"""

from __future__ import annotations

import argparse
import logging
import sys

import pandas as pd
import yaml

log = logging.getLogger("drought")

STAGE_COMMANDS = ("condition", "sensitivity", "anomaly", "cuts", "vhi")


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _years(text: str) -> list[int]:
    """"2019-2026" or "2024,2026" or "2026"."""
    out: list[int] = []
    for part in str(text).split(","):
        if "-" in part:
            a, b = part.split("-")
            out += list(range(int(a), int(b) + 1))
        elif part.strip():
            out.append(int(part))
    return sorted(set(out))


def _aois(a) -> list[str] | None:
    from drought.core import aoi
    return [aoi.resolve(x) for x in a.aoi] if getattr(a, "aoi", None) else None


def _one_aoi(a) -> str:
    sites = _aois(a)
    if not sites or len(sites) != 1:
        sys.exit("drought: give exactly one --aoi")
    return sites[0]


# --------------------------------------------------------------------------- #
# Commands. Each imports what it needs only when it runs, so overrides are
# in place first.
# --------------------------------------------------------------------------- #
def _run(a) -> None:
    from drought.pipeline import runner
    kw = {"years_": _years(a.years) if a.years else None, "workers": a.workers}
    sites = _aois(a)
    if a.dry_run:
        for name, site, will in runner.plan(a.stage, sites, force=a.force, **kw):
            print(f"  {'run ' if will else 'done'}  {name}{' ' + site if site else ''}")
        return
    runner.run(a.stage, sites, force=a.force, **kw)


def _fetch(a) -> None:
    from drought.core.config import years
    from drought.pipeline.stages import STAGES, sites_for
    from drought.sources import fetch
    sites = sites_for(STAGES["fetch"], _aois(a))
    failed = fetch.run(sites, _years(a.years) if a.years else years(), workers=a.workers,
                       vhi=not a.no_vhi, dry_run=a.dry_run)
    if failed:
        sys.exit(1)


def _stage(a) -> None:
    from drought.pipeline.stages import STAGES, sites_for
    st = STAGES[a.command]
    for site in sites_for(st, _aois(a)):
        log.info("== %s %s", a.command, site)
        st.run(site)


def _mask(a) -> None:
    from drought import mask
    from drought.core import config, results
    if a.action == "build":
        from drought.pipeline.stages import STAGES
        STAGES["mask"].run(None)
        return
    master = mask.master_parcels(a.shp_dir, verbose=False)
    from drought.viz import mask as mask_fig
    if a.action == "check":
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
    elif a.action == "figures":
        areas = {n: float(master[master["_code"].isin(c)].geometry.area.sum() / 1e4)
                 for n, c in mask.GROUPS.items()}
        mask_fig.figure_psf()
        mask_fig.figure_parcel_size(master=master)
        mask_fig.area_check(master=master)
        for site in _aois(a) or ["linthebene"]:
            mask_fig.figure_purity_explainer(site=site, areas=areas)
            mask_fig.figure_purity_map(site, group="intensive")
    else:                                               # overview
        from drought.sources import rgb
        sites = _aois(a) or list(config.CFG.condition.sites)
        for s in sites:
            if not rgb.path(s).exists():
                rgb.fetch(s)
        mask_fig.figure_mask_overview(sites=sites, thresholds=tuple(a.thresholds),
                                      master=master)


def _evaluate(a) -> None:
    from drought.pipeline.report import evaluate
    evaluate(_one_aoi(a))


def _precip(a) -> None:
    from drought.core.config import CFG, station_dir
    if a.download or not any(station_dir().glob("ogd-smn_*_d_*.csv")):
        from drought.sources import meteoswiss
        meteoswiss.download(station_dir(), CFG.precip.station_name)
    from drought.pipeline.report import report
    report(parts=("precip",))


def _report(a) -> None:
    from drought.pipeline.report import report
    parts = tuple(a.only) if a.only else ("condition", "sensitivity", "anomaly", "vhi",
                                          "precip")
    report(_aois(a), parts=parts)


def _inspect(a) -> None:
    import drought as dr
    from drought.analysis import detect, evaluate, labels, metrics
    site = _one_aoi(a)
    if a.action == "labels":
        lab = labels.load_labels(site)
        t = metrics.labels_last_cut_paired(lab, end=a.end)
        if t.empty:
            sys.exit("no parcel labelled in the current year and a baseline year")
        cur = int(dr.CFG.project.current_year)
        print(f"hand labels only, every year cut at DOY {a.end}:")
        print(f"  {len(t)} parcels labelled in {cur} and >= 1 baseline year "
              f"(median {t.n_baseline_years.median():.0f} baseline years each)")
        print(f"  baseline mean last cut  DOY {t.baseline.mean():.1f}")
        print(f"  {cur} mean last cut      DOY {t.current.mean():.1f}")
        print(f"  median shift            {t.shift_days.median():+.1f} days")
        print(f"  parcels stopping sooner {int((t.shift_days < 0).sum())} of {len(t)}")
        return
    df = dr.load_site(site, verbose=a.action == "detect", labelled=a.labelled)
    if a.action == "detect":
        print(metrics.summarise(detect.detect_all(df)).to_string())
    elif a.action == "errors":
        lab = labels.check_labels(labels.load_labels(site), df)
        by_year = evaluate.score_by_year(lab, df)
        cur = int(dr.CFG.project.current_year)
        rest = by_year.drop(index=cur, errors="ignore")[["tp", "fp", "fn"]].sum()
        p = rest.tp / (rest.tp + rest.fp)
        r = rest.tp / (rest.tp + rest.fn)
        print(f"\nwithout {cur}: precision {p:.3f}  recall {r:.3f}  F1 {2 * p * r / (p + r):.3f}")
        if cur in by_year.index:
            y = by_year.loc[cur]
            print(f"{cur} alone: {int(y.tp)} hits, {int(y.fp)} false positives, "
                  f"{int(y.fn)} misses")
        fp = evaluate.false_positive_positions(lab, df)
        if not fp.empty:
            fp["window"] = pd.cut(fp.doy, bins=list(range(60, 320, 20)), right=False)
            tab = fp.pivot_table(index="window", columns=fp.year == cur, values="doy",
                                 aggfunc="count", fill_value=0, observed=False)
            tab.columns = [str(cur) if c else "other years" for c in tab.columns]
            print("\nfalse positives by day of year:")
            print(tab.to_string())
        miss = evaluate.miss_positions(lab, df)
        if not miss.empty:
            print(f"\nmisses: {len(miss)}, of them last cuts of their season: "
                  f"{int(miss.is_last.sum())}")
    else:                                               # field
        from drought.viz import cuts as cut_fig
        if not a.parcel:
            lab = labels.load_labels(site)
            years = set(a.years)
            by_pid: dict[str, set] = {}
            for key in lab:
                pid, yr = key.rsplit("|", 1)
                by_pid.setdefault(pid, set()).add(int(yr))
            hits = sorted(pid for pid, ys in by_pid.items() if years <= ys)
            print(f"labelled in all of {sorted(years)}: {len(hits)} field(s)")
            for pid in hits:
                print(f"  {pid}")
            return
        ev = detect.detect_all(df, verbose=False)
        cut_fig.plot_parcel(df, a.parcel, ev, site=site, years=a.years, save=True)


def _demo(a) -> None:
    from drought import demo
    demo.run(fresh=a.fresh)


def _synthetic(a) -> None:
    from drought import synthetic
    synthetic.generate(seed=a.seed, scale=a.scale)


def _notebook(a) -> None:
    from drought import notebook
    notebook.run(check=a.check)


def _gallery(a) -> None:
    from drought.viz import gallery
    gallery.build()


def _provenance(a) -> None:
    import json

    from drought.core import provenance
    try:
        for f in a.files:
            if a.action == "verify":
                provenance.verify(f, a.kind)
                print(f"  ok  {f}")
            else:
                print(json.dumps(json.loads(provenance.sidecar_path(f).read_text()), indent=1))
    except provenance.ProvenanceError as exc:
        sys.exit(f"drought: {exc}")


def _results_show(a) -> None:
    import json

    from drought.core import results
    try:
        doc = results.load(strict=not a.allow_stale)
    except FileNotFoundError:
        sys.exit(f"no results yet: {results.path()}")
    except RuntimeError as exc:
        sys.exit(f"drought: {exc}")
    print(f"# {results.path()}")
    print(json.dumps({k: v["data"] for k, v in doc["sections"].items()}, indent=1))


def _config_show(a) -> None:
    from drought.core import config
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


def _aoi_arg(p, required=False) -> None:
    p.add_argument("--aoi", action="append", required=required, metavar="AOI",
                   help="site name from config.yaml or a GeoJSON file (repeatable)")


def build_parser() -> argparse.ArgumentParser:
    """The argument parser for the `drought` command."""
    from drought.core.provenance import KINDS
    from drought.pipeline.stages import STAGES

    p = argparse.ArgumentParser(
        prog="drought", description="Reading a drought from space: the pipeline.",
        epilog="Precedence: config.yaml < --config < --synthetic < --set.")
    p.add_argument("--config", metavar="FILE", help="use this YAML instead of config.yaml")
    p.add_argument("--synthetic", action="store_true",
                   help="use the synthetic data profile (data/synthetic, output/synthetic)")
    p.add_argument("--set", dest="overrides", metavar="KEY=VALUE", action="append",
                   type=_parse_set, default=[],
                   help="override one setting, e.g. --set mask.purity_min=0.8 (repeatable)")
    p.add_argument("-q", "--quiet", action="store_true", help="warnings and errors only")
    p.add_argument("-v", "--verbose", action="store_true", help="debug output")
    sub = p.add_subparsers(dest="command", required=True, metavar="COMMAND")

    c = sub.add_parser("run", help="a stage and everything it needs that is not done")
    c.add_argument("stage", nargs="?", default="report", choices=list(STAGES))
    _aoi_arg(c)
    c.add_argument("--years", help="for fetch: e.g. 2019-2026 (default: baseline + current)")
    c.add_argument("--workers", type=int, help="for fetch: parallel downloads")
    c.add_argument("--force", action="store_true", help="re-run the stage even if done")
    c.add_argument("--dry-run", action="store_true", help="list what would run")
    c.set_defaults(func=_run)

    c = sub.add_parser("fetch", help=STAGES["fetch"].help)
    _aoi_arg(c)
    c.add_argument("--years", help="e.g. 2019-2026 or 2024,2026 (default: baseline + current)")
    c.add_argument("--workers", type=int, help="parallel downloads (default swisseo.workers)")
    c.add_argument("--no-vhi", action="store_true", help="skip the VHI maps")
    c.add_argument("--dry-run", action="store_true", help="count what exists, download nothing")
    c.set_defaults(func=_fetch)

    m = sub.add_parser("mask", help="the grassland purity mask").add_subparsers(
        dest="action", required=True, metavar="ACTION")
    for name, text in (("build", STAGES["mask"].help),
                       ("check", "area cost per threshold; PSF model check"),
                       ("figures", "the mask figures"),
                       ("overview", "imagery, field outlines and purity, one row per site")):
        c = m.add_parser(name, help=text)
        c.add_argument("--shp-dir", help="folder with the yearly shapefiles (default: config)")
        if name == "check":
            c.add_argument("--site", help="restrict to one site window")
        if name in ("figures", "overview"):
            _aoi_arg(c)
        if name == "overview":
            c.add_argument("--thresholds", nargs="+", type=float, default=[0.0, 0.5, 0.9],
                           help="purity thresholds, one column each (default: 0 0.5 0.9)")
        c.set_defaults(func=_mask)

    for name in STAGE_COMMANDS:
        c = sub.add_parser(name, help=STAGES[name].help)
        _aoi_arg(c)
        c.set_defaults(func=_stage)

    c = sub.add_parser("evaluate", help=STAGES["evaluate"].help)
    _aoi_arg(c, required=True)
    c.set_defaults(func=_evaluate)

    c = sub.add_parser("precip", help="cumulative precipitation and P - ET0")
    c.add_argument("--download", action="store_true", help="fetch the station series again")
    c.set_defaults(func=_precip)

    c = sub.add_parser("report", help=STAGES["report"].help)
    _aoi_arg(c)
    c.add_argument("--only", action="append",
                   choices=["condition", "sensitivity", "anomaly", "vhi", "precip"])
    c.set_defaults(func=_report)

    ins = sub.add_parser("inspect", help="look closer at the cut detection").add_subparsers(
        dest="action", required=True, metavar="ACTION")
    for name, text in (("detect", "detected cuts, per-year summary"),
                       ("errors", "where the false positives and misses fall"),
                       ("field", "one field's series with its detected cuts (without "
                        "--parcel: list fields labelled in all --years)"),
                       ("labels", "last-cut shift from the hand labels alone")):
        c = ins.add_parser(name, help=text)
        _aoi_arg(c, required=True)
        if name in ("detect", "errors", "field"):
            c.add_argument("--labelled", action="store_true",
                           help="the labelled parcels' series instead of the sample")
        if name == "field":
            c.add_argument("--parcel", help="parcel_id; omit to list candidates")
            c.add_argument("--years", nargs="+", type=int, default=[2025, 2026])
        if name == "labels":
            c.add_argument("--end", type=int, default=228, help="comparison end DOY")
        c.set_defaults(func=_inspect)

    c = sub.add_parser("demo", help="the pipeline on synthetic data, self-checking")
    c.add_argument("--fresh", action="store_true", help="regenerate the synthetic inputs")
    c.set_defaults(func=_demo)

    syn = sub.add_parser("synthetic", help="synthetic inputs").add_subparsers(
        dest="action", required=True, metavar="ACTION")
    c = syn.add_parser("generate", help="write every synthetic input")
    c.add_argument("--seed", type=int, help="master seed (default: config)")
    c.add_argument("--scale", type=float, help="share of the real parcel count")
    c.set_defaults(func=_synthetic)
    c = sub.add_parser("notebook", help="run notebooks/walkthrough.ipynb on synthetic data")
    c.add_argument("--check", action="store_true", help="run a copy; fail on any error (CI)")
    c.set_defaults(func=_notebook)
    sub.add_parser("gallery", help="refresh docs/figures/model").set_defaults(func=_gallery)

    pv = sub.add_parser("provenance", help="provenance records of data files").add_subparsers(
        dest="action", required=True, metavar="ACTION")
    for name, text in (("verify", "check files against the current configuration and mask"),
                       ("show", "print a file's record")):
        c = pv.add_parser(name, help=text)
        c.add_argument("files", nargs="+")
        if name == "verify":
            c.add_argument("--kind", required=True, choices=list(KINDS))
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
    from drought.core.log import setup
    setup(logging.WARNING if args.quiet else logging.DEBUG if args.verbose else logging.INFO)
    import matplotlib
    matplotlib.use("Agg")                         # a command line never opens windows

    from drought.core import config
    synthetic = args.synthetic or args.command == "demo"
    try:
        if args.config:
            config.use_config(args.config)
        if synthetic:
            config.activate("synthetic")
        for key, value in args.overrides:
            config.override(key, value)
    except (KeyError, ValueError) as exc:
        sys.exit(f"drought: {exc.args[0]}")
    if synthetic:
        from drought.synthetic import stamp_figures
        stamp_figures()
    from drought.core.provenance import ProvenanceError
    try:
        args.func(args)
    except ProvenanceError as exc:
        sys.exit(f"drought: refused — {exc}")
    except FileNotFoundError as exc:
        sys.exit(f"drought: {exc.args[0] if exc.args else exc}")


if __name__ == "__main__":
    main()
