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
    from . import mask
    mask.build_mask(shp_dir=a.shp_dir)


def _mask_check(a) -> None:
    from . import mask
    master = mask.master_parcels(a.shp_dir, verbose=False)
    mask.area_table(master, site=a.site)
    mask.psf_model_check(master)


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
    from . import trajectory
    res = trajectory.plot(site=a.site, path=a.path, index=a.index)
    if res is not None:
        print(res[["site", "regime", "integral_pct", "rank", "upto", "n_days"]]
              .to_string(index=False))


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
    label = "Synthetic station" if config.active_profile() == "synthetic" else None
    precip.main(daily=daily, abbr=abbr, lat=lat, station_label=label)


def _cuts_detect(a) -> None:
    import drought as dr
    df = dr.load_site(a.site)
    ev = dr.detect.detect_all(df)
    print(dr.metrics.summarise(ev).to_string())


def _cuts_score(a) -> None:
    import drought as dr

    from . import changepoint, evaluate, labels
    df = dr.load_site(a.site, verbose=False)
    lab = labels.check_labels(labels.load_labels(a.site), df)
    evaluate.assert_wired()
    print("heuristic (frozen configuration):")
    evaluate.score(lab, df)
    print(evaluate.score_by_year(lab, df)[["n_parcel_years", "precision", "recall"]]
          .to_string())
    cp = changepoint.score_detector(lab, df, changepoint.detect_changepoints)
    print(f"\nchangepoint alternative: precision {cp['precision']}  recall "
          f"{cp['recall']}  F1 {cp['f1']}")


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
                             ("figures", _mask_figures, "the mask figures")):
        c = m.add_parser(name, help=text)
        c.add_argument("--shp-dir", help="folder with the yearly shapefiles (default: config)")
        if name == "check":
            c.add_argument("--site", help="restrict to one site window")
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

    c = sub.add_parser("precip", help="cumulative precipitation and P - ET0")
    c.add_argument("--download", action="store_true",
                   help="fetch from MeteoSwiss instead of reading paths.station")
    c.set_defaults(func=_precip)

    cu = sub.add_parser("cuts", help="cut detection").add_subparsers(
        dest="action", required=True, metavar="ACTION")
    for name, func, text in (("detect", _cuts_detect, "detect cuts, per-year summary"),
                             ("score", _cuts_score, "score against the hand labels")):
        c = cu.add_parser(name, help=text)
        c.add_argument("--site", required=True)
        c.set_defaults(func=func)

    sub.add_parser("gallery", help="refresh docs/figures/{model,demo}").set_defaults(func=_gallery)

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
    args.func(args)


if __name__ == "__main__":
    main()
