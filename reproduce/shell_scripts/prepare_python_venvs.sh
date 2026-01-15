#!/usr/bin/bash
set -euo pipefail

readarray -t venvs < <(
	yq -r '.paths.core.venv.py_models_dir' "$1"
	yq -r '.paths.core.venv.analysis_dir' "$1"
)

for path in "${venvs[@]}"; do
	python3 -m venv "$path"
	source "${path}bin/activate"
	python3 -m pip install -r "${path}requirements.txt"
done
