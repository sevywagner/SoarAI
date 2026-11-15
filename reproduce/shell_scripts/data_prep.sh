#!/usr/bin/bash
set -euo pipefail

../build/src/data_prep $1

run_root=`yq -r '.paths.core.run_root' $1`
run_id=`yq -r '.run.run_id' $1`

cp $1 "${run_root}${run_id}/"
