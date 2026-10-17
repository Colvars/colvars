#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="${OPENMM_COLVARS_BUILD_DIR:-${repo_root}/build-openmm}"
install_dir="${OPENMM_COLVARS_INSTALL_DIR:-${repo_root}/install-openmm}"

if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "No active conda or pixi environment with OpenMM detected." >&2
    echo "Run this script from the OpenMM pixi environment or activate a compatible OpenMM conda env first." >&2
    exit 1
fi

openmm_source_dir="${OPENMM_SOURCE:-}"
if [[ -z "${openmm_source_dir}" && -d "${repo_root}/openmm-source" ]]; then
    openmm_source_dir="${repo_root}/openmm-source"
fi

lepton_dir="${LEPTON_DIR:-}"
if [[ -z "${lepton_dir}" && -n "${openmm_source_dir}" && -d "${openmm_source_dir}/libraries/lepton" ]]; then
    lepton_dir="${openmm_source_dir}/libraries/lepton"
fi

cmake_args=(
    -S "${repo_root}/openmm"
    -B "${build_dir}"
    -DCMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE:-Release}"
    -DCMAKE_INSTALL_PREFIX="${install_dir}"
    -DCMAKE_PREFIX_PATH="${CONDA_PREFIX}"
    -DCOLVARS_SOURCE_DIR="${repo_root}"
    -DOPENMM_COLVARS_BUILD_PYTHON=ON
    -DOPENMM_COLVARS_ENABLE_CUDA_KERNELS="${OPENMM_COLVARS_ENABLE_CUDA_KERNELS:-OFF}"
)

if [[ -n "${openmm_source_dir}" ]]; then
    cmake_args+=( -DOPENMM_SOURCE_DIR="${openmm_source_dir}" )
fi
if [[ -n "${lepton_dir}" && -d "${lepton_dir}" ]]; then
    cmake_args+=( -DLEPTON_DIR="${lepton_dir}" )
fi

mkdir -p "${build_dir}"
cmake "${cmake_args[@]}"
cmake --build "${build_dir}" --parallel "$(nproc --all)"
cmake --install "${build_dir}"

python - <<'PY'
import openmm
from openmm import Platform, System
from openmm_colvars import ColvarsForce
from openmm_fep import FEPAtoms, FEPConfig

print("OpenMM:", openmm.__version__)
print("OpenMM platforms:", [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())])

system = System()
system.addParticle(1.0)
system.addParticle(1.0)
force = ColvarsForce.fromConfigString("colvarsTrajFrequency 0\n", "colvars_ci")
system.addForce(force)
print("ColvarsForce:", type(force).__name__)
print("FEP API:", FEPAtoms, FEPConfig)
PY
