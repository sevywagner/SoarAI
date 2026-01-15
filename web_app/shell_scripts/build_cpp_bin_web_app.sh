#!/usr/bin/bash
set -euo pipefail

mkdir -p -- ../build
cd ../build
cmake -DSOAR_BUILD_API=ON ..
make
cd ../web_app
