# Fresh clone:   make install && make demo
PY ?= python

.PHONY: install demo demo-fresh gallery notebook site test lint check clean

install:            ## editable install with test and lint tools
	$(PY) -m pip install -e ".[dev]"

demo:               ## whole pipeline on synthetic data -> output/synthetic/
	PYTHONPATH=src $(PY) -m drought.cli demo

demo-fresh:         ## same, regenerating the synthetic inputs first
	PYTHONPATH=src $(PY) -m drought.cli demo --fresh

gallery:            ## refresh docs/figures/model (the sensor-response figure)
	PYTHONPATH=src $(PY) -m drought.cli gallery

notebook:           ## re-run notebooks/walkthrough.ipynb on synthetic data (commit the result)
	PYTHONPATH=src $(PY) -m drought.cli notebook

site:               ## build the website into site/_build/html (pip install "jupyter-book>=1,<2" sphinxcontrib-mermaid)
	rm -rf site/src && mkdir -p site/src
	cp -r README.md LICENSE docs notebooks site/src/
	cp site/_config.yml site/_toc.yml site/src/
	jupyter-book build site/src --path-output site

test:
	$(PY) -m pytest

lint:
	$(PY) -m ruff check src tests

check: lint test demo

clean:              ## remove synthetic inputs and outputs (never touches real data)
	rm -rf data/synthetic output/synthetic
