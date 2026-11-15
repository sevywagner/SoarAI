#!/usr/bin/bash
set -euo pipefail

config=../configs/soar.yaml

./shell_scripts/validate_config.sh $config
./shell_scripts/download_h5.sh $config
./shell_scripts/prepare_python_venvs.sh $config
./shell_scripts/build_cpp_bin.sh
./shell_scripts/data_prep.sh $config
echo "========================================"
echo "Training XGBoost"
echo "========================================"
./shell_scripts/train.sh $config
echo "========================================"
echo "Training Neural Network"
echo "========================================"
./shell_scripts/create_py_models.sh $config
echo "========================================"
echo "Analysis"
echo "========================================"
./shell_scripts/analysis.sh $config
