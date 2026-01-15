#!/usr/bin/bash
set -euo pipefail

path=`yq -r '.paths.core.venv.analysis_dir' ../configs/soar.yaml`

source "${path}bin/activate"
python3 "${path}analysis.py" ../configs/soar.yaml NH4
python3 "${path}analysis.py" ../configs/soar.yaml NO
