"""The committed walkthrough must come from a synthetic run and leak no local path."""

import json
import re

from drought.notebook import NOTEBOOK

NB = json.loads(NOTEBOOK.read_text())


def _outputs():
    for c in NB["cells"]:
        for o in c.get("outputs", []):
            yield "".join(o.get("text", "")) + "".join(o.get("data", {}).get("text/plain", ""))


def test_walkthrough_is_set_to_synthetic():
    first = next("".join(c["source"]) for c in NB["cells"] if c["cell_type"] == "code")
    assert re.search(r"^REAL\s*=\s*False", first, re.M)


def test_committed_outputs_are_from_a_synthetic_run():
    text = "\n".join(_outputs())
    assert "profile: synthetic" in text, "run `make notebook` before committing"
    assert "profile: real" not in text, "a notebook run on real data must never be committed"


def test_no_local_paths_in_outputs():
    for t in _outputs():
        assert not re.search(r"/Users/|/home/|[A-Z]:\\\\Users", t), t[:200]


def test_every_code_cell_ran_without_error():
    for c in NB["cells"]:
        if c["cell_type"] == "code":
            assert c.get("execution_count"), "unexecuted cell: run `make notebook`"
            assert not any(o["output_type"] == "error" for o in c.get("outputs", []))
