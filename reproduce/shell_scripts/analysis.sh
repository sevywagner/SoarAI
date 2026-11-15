#!/usr/bin/bash
set -euo pipefail

path=`yq -r '.paths.core.venv.analysis_dir' $1`
run_root=`yq -r '.paths.core.run_root' $1`
run_dir=${run_root}`yq -r '.run.run_id' $1`/
data_dir=${run_dir}`yq -r '.paths.core.data_dir' $1`
h5_file=`yq -r '.paths.analysis.data_spectra_h5' $1`

cp "../data/${h5_file}" "${data_dir}/${h5_file}"

source "${path}bin/activate"
python3 "${path}analysis.py" $1 final
python3 "${path}van_krevelen.py" $1
python3 "${path}mass_defect.py" $1
python3 "${path}plot_nh4_peak_fit.py" $1
