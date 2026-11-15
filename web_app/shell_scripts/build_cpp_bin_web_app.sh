#!/usr/bin/bash
set -euo pipefail

mkdir -p -- ../build
cd ../build
cmake -DSOAR_BUILD_API=ON ..
make
cd ../web_app

model=$(yq -r '.paths.core.run_root + .run.run_id + "/" + .paths.core.model_output_dir + "xgboost.cmod"' "$1")
api_model=$(yq -r '.paths.api.model_path' "$1")

[[ "$model" = /* ]] || model="../reproduce/$model"
mkdir -p -- "$(dirname -- "$api_model")"

if [[ -f "$model" ]]; then
    cp -- "$model" "$api_model"
else
    echo "Please place a trained XGBoost model at: $api_model"
fi

[[ -f "$api_model" ]] || exit 1
