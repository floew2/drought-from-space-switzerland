# Fresh clone:   make install && make demo
PY ?= python

.PHONY: install demo demo-fresh gallery test lint check clean

install:            ## editable install with test and lint tools
	$(PY) -m pip install -e ".[dev]"

demo:               ## whole pipeline on synthetic data -> output/synthetic/
	PYTHONPATH=src $(PY) -m drought.cli demo

demo-fresh:         ## same, regenerating the synthetic inputs first
	PYTHONPATH=src $(PY) -m drought.cli demo --fresh

gallery: demo       ## refresh the committed figures in docs/figures/{model,demo}
	PYTHONPATH=src $(PY) -m drought.cli gallery

test:
	$(PY) -m pytest

lint:
	$(PY) -m ruff check src tests

check: lint test demo

clean:              ## remove synthetic inputs and outputs (never touches real data)
	rm -rf data/synthetic output/synthetic
