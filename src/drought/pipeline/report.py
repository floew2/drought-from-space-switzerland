"""The last two stages: scoring the cut detector, and the report.

    evaluate(site)    precision, recall and F1 against the hand labels
                      -> results.json `detector.<site>`
    report(sites)     every figure, and the headline numbers -> results.json:
                      condition, sensitivity, anomaly, vhi, precipitation

Each part runs if its inputs exist and says so if they do not, so a report on
partial data is still a report.
"""

from __future__ import annotations

import logging

from drought.core import results
from drought.core.config import CFG, active_profile, series_dir, station_dir

__all__ = ["evaluate", "report"]

log = logging.getLogger(__name__)


def evaluate(site: str, labelled: bool = True) -> dict:
    """Score the frozen detector and the changepoint alternative on one site."""
    import drought as dr
    from drought.analysis import changepoint, labels
    from drought.analysis import evaluate as ev

    df = dr.load_site(site, verbose=False, labelled=labelled)
    lab = labels.check_labels(labels.load_labels(site), df)
    ev.assert_wired()
    log.info("heuristic (frozen configuration):")
    h = ev.score(lab, df)
    by_year = ev.score_by_year(lab, df)
    log.info("%s", by_year[["n_parcel_years", "precision", "recall"]].to_string())
    cp = changepoint.score_detector(lab, df, changepoint.detect_changepoints)
    log.info("\nchangepoint alternative: precision %s  recall %s  F1 %s",
             cp["precision"], cp["recall"], cp["f1"])
    out = {"labelled_parcel_years": len(lab), "precision": h["precision"],
           "recall": h["recall"], "f1": h["f1"],
           "recall_per_year": [float(by_year.recall.min()), float(by_year.recall.max())],
           "changepoint_f1": cp["f1"]}
    pattern = "cuts_*_labelled.csv" if labelled else "cuts_????.csv"
    results.update(f"detector.{site}", out, inputs=sorted(series_dir(site).glob(pattern)))
    return out


def _condition(sites) -> None:
    from drought.viz import condition
    res = condition.plot(sites)
    if res is None:
        return
    log.info("%s", res[["site", "regime", "integral_pct", "rank", "upto", "n_days"]]
             .to_string(index=False))
    results.update("condition", {
        f"{r.site}/{r.regime}": {"index": CFG.detect.index,
                                 "integral_pct": round(float(r.integral_pct), 1),
                                 "rank": int(r.rank), "n_baseline_years": int(r.n_baseline_years),
                                 "upto_doy": int(r.upto), "shared_days": int(r.n_days)}
        for r in res.itertuples()},
        inputs=[series_dir(s) / "condition.csv" for s in res.site.unique()])


def _sensitivity(sites) -> None:
    from drought.viz import sensitivity
    out, used = {}, CFG.mask.purity_min
    for site in sites:
        for regime in CFG.analysis.regimes:
            log.info("\n%s / %s", site, regime)
            t = sensitivity.figure(site, regime)
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
        log.info("\n  spread of the anomaly over thresholds 0.5-0.95, as %% of the anomaly:\n%s",
                 "\n".join(f"    {k:<24}{v['spread_0p5_up_pct_of_anomaly']:>6.1f} %"
                           for k, v in out.items()))


def _anomaly(sites) -> None:
    from drought.viz import anomaly
    out = {}
    for site in sites:
        r = anomaly.figure(site)
        r.pop("file", None)
        out[site] = r
    if out:
        results.update("anomaly", out, inputs=[series_dir(s) / "anomaly.csv" for s in out])


def _vhi(sites) -> None:
    from drought.viz import vhi
    out = vhi.figure(sites)
    if out:
        results.update("vhi", out, inputs=[series_dir(s) / "vhi_series.csv" for s in sites])


def _precip() -> None:
    from drought.sources import meteoswiss
    from drought.viz import precip
    daily, abbr, lat = meteoswiss.load_local(station_dir())
    label = "Synthetic station" if active_profile() == "synthetic" else None
    res = precip.figure(daily, abbr, lat, station_label=label)
    res.pop("file", None)
    results.update("precipitation", res,
                   inputs=sorted(station_dir().glob("ogd-smn_*_d_*.csv")))


def report(sites: list[str] | None = None, parts=("condition", "sensitivity", "anomaly",
                                                   "vhi", "precip")) -> list[str]:
    """Draw every figure and write the headline numbers. Returns the parts skipped."""
    pick = {"condition": lambda: _condition(sites or list(CFG.condition.sites)),
            "sensitivity": lambda: _sensitivity(sites or list(CFG.condition.sites)),
            "anomaly": lambda: _anomaly(sites or list(CFG.anomaly.sites)),
            "vhi": lambda: _vhi(sites or list(CFG.condition.sites)),
            "precip": _precip}
    skipped = []
    for name in parts:
        log.info("\n=== report: %s %s", name, "=" * (50 - len(name)))
        try:
            pick[name]()
        except FileNotFoundError as exc:
            log.warning("  [skipped] %s", exc)
            skipped.append(name)
    return skipped
