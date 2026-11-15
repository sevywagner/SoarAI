#!/usr/bin/bash
set -eou pipefail

run_root=`yq -r '.paths.core.run_root' ../configs/soar.yaml`
run_id=`yq -r '.run.run_id' ../configs/soar.yaml`
run_dir="${run_root}${run_id}/"

./shell_scripts/download_h5.sh ../configs/soar.yaml

mkdir -p -- $run_dir
sudo docker build -t "soarai:1.1.0" ../ -f ./docker/dockerfile
sudo docker run --mount type=bind,source="${run_dir}",target="/app/reproduce/${run_dir}" soarai:1.1.0
