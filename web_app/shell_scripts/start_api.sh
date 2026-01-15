#!/usr/bin/bash
set -euo pipefail

./shell_scripts/build_cpp_bin_web_app.sh

path=`yq -r '.paths.api.peak_fit_binary' ../configs/soar.yaml`

if [[ -f $path ]]; then
	echo "$path exists. Skipping build. "
else
	./shell_scripts/build_peak_fit_bin.sh
fi

../build/src/api_driver ../configs/soar.yaml NH4 &
../build/src/api_driver ../configs/soar.yaml NO &
