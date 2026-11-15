#!/usr/bin/bash
set -euo pipefail

path=`yq -r '.paths.core.venv.analysis_dir' $1`
run_root=`yq -r '.paths.core.run_root' $1`
run_dir=${run_root}`yq -r '.run.run_id' $1`/
data_dir=${run_dir}`yq -r '.paths.core.data_dir' $1`
h5_file=`yq -r '.paths.analysis.data_spectra_h5' $1`

h5_path="../data/$h5_file"

if [[ -f $h5_path ]]; then
    echo "$path exists. Skipping h5 download. "
else
    cd ../data/
    wget https://zenodo.org/records/24214638/files/a.h5
fi
