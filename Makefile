.PHONY: run test clean

run:
	./run_pipeline.sh -o results

test:
	python3 -m pytest -q tests

clean:
	rm -rf results .pytest_cache tests/__pycache__ scripts/__pycache__
