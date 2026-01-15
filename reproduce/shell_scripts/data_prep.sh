#!/usr/bin/bash
set -euo pipefail

../build/src/data_prep ../configs/soar.yaml

run_root=`yq -r '.paths.core.run_root' ../configs/soar.yaml`
run_id=`yq -r '.run.run_id' ../configs/soar.yaml`

cp ../configs/soar.yaml "${run_root}${run_id}/"
