#!/usr/bin/bash
set -euo pipefail

../build/src/train ../configs/soar.yaml NH4
../build/src/train ../configs/soar.yaml NO
