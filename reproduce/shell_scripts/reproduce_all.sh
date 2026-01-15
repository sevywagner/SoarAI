#!/usr/bin/bash
set -euo pipefail

./shell_scripts/validate_config.sh ../configs/soar.yaml -a
./shell_scripts/prepare_python_venvs.sh ../configs/soar.yaml
./shell_scripts/data_prep.sh
./shell_scripts/train.sh
./shell_scripts/create_py_models.sh
./shell_scripts/analysis.sh
