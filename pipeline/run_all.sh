#!/usr/bin/env bash
# Re-run network-dependent inputs and the public downstream analysis.
# Archived external B-cell predictor outputs are not recreated by this script.
set -euo pipefail
cd "$(dirname "$0")/.."
PY=${PY:-python3}

echo "== 1/9 conservation sequence panel =="
$PY pipeline/fetch_conservation_sequences.py

echo "== 2/9 population HLA panel =="
$PY pipeline/download_population_panel.py

echo "== 3/9 class-I predictions =="
$PY pipeline/run_predictions.py --class 1

echo "== 4/9 class-II predictions =="
$PY pipeline/run_predictions.py --class 2

echo "== 5/9 conservation analysis =="
$PY pipeline/conservation.py

echo "== 6/9 human-proteome sequence similarity =="
$PY pipeline/pepmatch_screen.py

echo "== 7/9 core collapse, scoring, and selected set =="
$PY pipeline/collapse_and_score.py

echo "== 8/9 genotype coverage and B-cell table post-processing =="
$PY pipeline/population_coverage.py
$PY pipeline/bcell_processing.py

echo "== 9/9 master file and figures =="
$PY pipeline/build_masterfile.py
$PY paper/scripts/make_figures.py
$PY tools/build_manifest.py

echo "Done. Run 'make validate' for static release-integrity checks."
