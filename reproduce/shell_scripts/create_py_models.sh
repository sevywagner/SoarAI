#!/usr/bin/bash
set -euo pipefail

path=`yq -r '.paths.core.venv.py_models_dir' $1`

source "${path}bin/activate"
python3 "${path}src/create_models.py" -c $1 -n -d final
