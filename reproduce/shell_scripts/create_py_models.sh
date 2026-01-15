#!/usr/bin/bash
set -euo pipefail

path=`yq -r '.paths.core.venv.py_models_dir' ../configs/soar.yaml`

source "${path}bin/activate"
python3 "${path}src/create_models.py" -c ../configs/soar.yaml -r NH4 -nlg
python3 "${path}src/create_models.py" -c ../configs/soar.yaml -r NO -nlg
