"""The `drought` command: precedence, typo protection, and overrides that reach the code.

Each case runs in a fresh interpreter, because some modules copy settings into
constants at import and the test session has already imported them.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def _run(*args: str, code: str | None = None) -> subprocess.CompletedProcess:
    cmd = [sys.executable, "-c", code] if code else [sys.executable, "-m", "drought.cli", *args]
    env = {**os.environ, "PYTHONPATH": str(ROOT / "src")}
    env.pop("DROUGHT_PROFILE", None)
    return subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT, env=env)


def _shown(*args: str) -> str:
    """Last line printed by `drought <args>`."""
    r = _run(*args)
    assert r.returncode == 0, r.stderr
    return r.stdout.splitlines()[-1]


def test_help_lists_every_command():
    from drought.cli import build_parser

    text = build_parser().format_help()
    for cmd in ("demo", "run-all", "synthetic", "mask", "trajectory", "precip", "cuts",
                "gallery", "config"):
        assert cmd in text


def test_precedence_config_file_then_set(tmp_path):
    cfg = yaml.safe_load((ROOT / "config.yaml").read_text())
    cfg["mask"]["purity_min"] = 0.7
    alt = tmp_path / "alt.yaml"
    alt.write_text(yaml.safe_dump(cfg))

    assert _shown("config", "show", "mask.purity_min") == "0.9"
    assert _shown("--config", str(alt), "config", "show", "mask.purity_min") == "0.7"
    assert _shown("--config", str(alt), "--set", "mask.purity_min=0.6",
                  "config", "show", "mask.purity_min") == "0.6"


def test_synthetic_profile_switches_paths():
    assert _shown("--synthetic", "config", "show", "paths.figures") == "output/synthetic/figures"


def test_unknown_key_is_refused_not_added():
    r = _run("--set", "mask.purity_mn=0.8", "config", "show")
    assert r.returncode != 0
    assert "no setting 'mask.purity_mn'" in r.stderr
    assert "purity_min" in r.stderr                    # names what does exist


def test_override_reaches_the_module_that_reads_it():
    r = _run(code="from drought import config; config.override('mask.purity_min', 0.8); "
                  "from drought.mask import settings; print(settings.PURITY_MIN)")
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip() == "0.8"


def test_override_after_import_is_refused():
    """Changing a setting a module has already copied would change nothing; refuse it."""
    r = _run(code="from drought.mask import settings; from drought import config; "
                  "config.override('mask.purity_min', 0.8)")
    assert r.returncode != 0
    assert "already read the configuration" in r.stderr


def test_overriding_the_frozen_detector_warns():
    r = _run("--set", "detect.min_drop_abs=0.1", "config", "show", "detect.min_drop_abs")
    assert r.returncode == 0
    assert r.stdout.splitlines()[-1] == "0.1"
    assert "frozen" in r.stderr
