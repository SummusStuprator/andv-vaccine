.PHONY: validate figures master manifest release-check full

validate:
	python3 tools/validate_release.py

figures:
	python3 paper/scripts/make_figures.py

master:
	python3 pipeline/build_masterfile.py

manifest:
	python3 tools/build_manifest.py

release-check: manifest validate

full:
	bash pipeline/run_all.sh
