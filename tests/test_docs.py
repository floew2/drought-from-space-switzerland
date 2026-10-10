"""Every image and relative link in the docs must resolve.

This is the pre-submission gate for docs/figures/real/: the README and gallery
reference the real-run figures, and this test stays red until they are present.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DOCS = [ROOT / "README.md", *sorted((ROOT / "docs").glob("*.md"))]
LINK = re.compile(r"!?\[[^\]]*\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)|<img[^>]+src=\"([^\"]+)\"")


def _targets(md: Path):
    for m in LINK.finditer(md.read_text(encoding="utf-8")):
        t = m.group(1) or m.group(2)
        if t.startswith(("http://", "https://", "mailto:", "#")):
            continue
        yield t.split("#")[0]


@pytest.mark.parametrize("md", DOCS, ids=lambda p: p.name)
def test_relative_links_resolve(md):
    missing = [t for t in _targets(md) if t and not (md.parent / t).exists()]
    assert not missing, f"{md.name} references missing files: {missing}"


@pytest.mark.parametrize("md", DOCS, ids=lambda p: p.name)
def test_no_unfilled_placeholders(md):
    text = md.read_text(encoding="utf-8")
    assert "{{" not in text, f"{md.name} still has a {{{{...}}}} placeholder"


def _slug(heading: str) -> str:
    """GitHub's anchor for a heading."""
    h = re.sub(r"[^\w\- ]", "", heading.strip().lower())
    return h.replace(" ", "-")


@pytest.mark.parametrize("md", DOCS, ids=lambda p: p.name)
def test_section_anchors_exist(md):
    text = md.read_text(encoding="utf-8")
    bad = []
    for target, anchor in re.findall(r"\]\(([^)#\s]*)#([\w\-]+)\)", text):
        doc = (md.parent / target) if target else md
        heads = {_slug(h) for h in re.findall(r"^#+ (.+)$", doc.read_text(encoding="utf-8"), re.M)}
        if anchor not in heads:
            bad.append(f"{target}#{anchor}")
    assert not bad, f"{md.name} links to missing sections: {bad}"
