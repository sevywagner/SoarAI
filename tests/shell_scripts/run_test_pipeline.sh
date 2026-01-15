#!/usr/bin/bash
set -euo pipefail

../reproduce/shell_scripts/validate_config.sh ../configs/test_config.yaml
../reproduce/shell_scripts/build_cpp_bin.sh
../reproduce/shell_scripts/prepare_python_venvs.sh ../configs/test_config.yaml
../build/src/data_prep ../configs/test_config.yaml
../build/src/train ../configs/test_config.yaml NH4 1
../build/src/train ../configs/test_config.yaml NH4 2

py_models_path=`yq -r '.paths.core.venv.py_models_dir' ../configs/test_config.yaml`

source "${py_models_path}bin/activate"
python3 "${py_models_path}src/create_models.py" -c ../configs/test_config.yaml -r NH4 -nlg -d 1
python3 "${py_models_path}src/create_models.py" -c ../configs/test_config.yaml -r NH4 -nlg -d 2

run_root=`yq -r '.paths.core.run_root' ../configs/test_config.yaml`
run_id=`yq -r '.run.run_id' ../configs/test_config.yaml`
pred_dir=`yq -r '.paths.core.pred_output_dir' ../configs/test_config.yaml`
logit_dir=`yq -r '.paths.core.pred_logit_output_dir' ../configs/test_config.yaml`
full_logit_dir="${run_root}${run_id}/${pred_dir}${logit_dir}"

models=("xgboost" "nn" "lgbm" "gbdt")
failed=0

echo "[*] Comparing logits in dir \"$full_logit_dir\""
for model in "${models[@]}"; do
	set +e
	diff -q "${full_logit_dir}NH4_${model}_logits_1.txt" "${full_logit_dir}NH4_${model}_logits_2.txt" > /dev/null
	diff_status=$?
	set -e
	if [ $diff_status -eq 0 ]; then
    		echo "[*] ${model} passed deterministic training test"
	else
		if [ $diff_status -eq 2 ]; then
			echo "diff error in ${model} logit comparison"
			exit 2
		else
			echo "Non determinism detected in ${model}"
			failed=1
		fi
	fi
done

cp ../configs/test_config.yaml "${run_root}${run_id}/"

exit $failed
