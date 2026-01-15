#!/usr/bin/bash
set -euo pipefail

path=`yq -r '.paths.api.venv.utils_dir' ../configs/soar.yaml`

python3 -m venv "$path"
source "${path}bin/activate"
python3 -m pip install -r "${path}requirements.txt"
pyinstaller --distpath "${path}dist" --workpath $path --specpath "${path}fit_peaks.spec" "${path}fit_peaks.py"
