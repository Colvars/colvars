#!/usr/bin/env bash
set -euo pipefail
source_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
build_dir=${OPENMM_COLVARS_BUILD_DIR:-"$source_dir/../build-openmm"}
python_exe=${PYTHON:-python}
cmake -S "$source_dir" -B "$build_dir" \
  -DCMAKE_BUILD_TYPE=Release \
  -DPython_EXECUTABLE="$(command -v "$python_exe")" \
  -DOPENMM_COLVARS_ENABLE_CUDA="${OPENMM_COLVARS_ENABLE_CUDA:-OFF}" "$@"
cmake --build "$build_dir" --parallel "${BUILD_JOBS:-2}"
cmake --install "$build_dir"
