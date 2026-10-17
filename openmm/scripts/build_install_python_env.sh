#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="${OPENMM_COLVARS_BUILD_DIR:-${repo_root}/build}"

if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "CONDA_PREFIX is not set. Activate your micromamba/conda OpenMM environment first." >&2
    exit 1
fi

colvars_source_dir="${COLVARS_SOURCE_DIR:-${repo_root}/tests_dependence_src/colvars}"
openmm_source_dir="${OPENMM_SOURCE_DIR:-${repo_root}/tests_dependence_src/openmm}"
lepton_dir="${LEPTON_DIR:-${openmm_source_dir}/libraries/lepton}"
enable_cuda="${OPENMM_COLVARS_ENABLE_CUDA_KERNELS:-ON}"

if [[ ! -f "${colvars_source_dir}/src/colvarmodule.cpp" ]]; then
    echo "Colvars source not found at ${colvars_source_dir}." >&2
    echo "Set COLVARS_SOURCE_DIR to a Colvars checkout." >&2
    exit 1
fi

if [[ "${enable_cuda}" == "ON" || "${enable_cuda}" == "1" || "${enable_cuda}" == "TRUE" ]]; then
    if [[ ! -f "${openmm_source_dir}/platforms/cuda/include/CudaContext.h" ]]; then
        echo "OpenMM CUDA internal headers not found at ${openmm_source_dir}." >&2
        echo "Set OPENMM_SOURCE_DIR to an OpenMM source tree, or set OPENMM_COLVARS_ENABLE_CUDA_KERNELS=OFF for CPU/Reference builds." >&2
        exit 1
    fi

    if [[ -z "${OPENMM_CUDA_LIBRARY:-}" ]]; then
        openmm_cuda_library=""
        for candidate in \
            "${CONDA_PREFIX}/lib/plugins/libOpenMMCUDA.so" \
            "${CONDA_PREFIX}/lib/libOpenMMCUDA.so" \
            "${CONDA_PREFIX}/Library/lib/plugins/OpenMMCUDA.lib" \
            "${CONDA_PREFIX}/Library/lib/OpenMMCUDA.lib"; do
            if [[ -f "${candidate}" ]]; then
                openmm_cuda_library="${candidate}"
                break
            fi
        done
        if [[ -z "${openmm_cuda_library}" ]]; then
            echo "OpenMM CUDA library was not found under CONDA_PREFIX=${CONDA_PREFIX}." >&2
            echo "Set OPENMM_CUDA_LIBRARY to libOpenMMCUDA.so/OpenMMCUDA.lib." >&2
            exit 1
        fi
    else
        openmm_cuda_library="${OPENMM_CUDA_LIBRARY}"
    fi

    if [[ -z "${OPENMM_COLVARS_CUDA_TOOLKIT_INCLUDE_DIR:-}" ]]; then
        cuda_include_dir=""
        for candidate in \
            "${CONDA_PREFIX}/include" \
            "${CONDA_PREFIX}/targets/x86_64-linux/include" \
            "${CUDA_HOME:-}/include" \
            "${CUDA_PATH:-}/include" \
            "/usr/local/cuda/include"; do
            if [[ -n "${candidate}" && -f "${candidate}/cuda.h" && -f "${candidate}/builtin_types.h" ]]; then
                cuda_include_dir="${candidate}"
                break
            fi
        done
        if [[ -z "${cuda_include_dir}" ]]; then
            echo "CUDA driver headers cuda.h and builtin_types.h were not found." >&2
            echo "Install CUDA development headers in the active environment, or set OPENMM_COLVARS_CUDA_TOOLKIT_INCLUDE_DIR." >&2
            exit 1
        fi
    else
        cuda_include_dir="${OPENMM_COLVARS_CUDA_TOOLKIT_INCLUDE_DIR}"
    fi
fi

cmake_args=(
    -S "${repo_root}"
    -B "${build_dir}"
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_PREFIX_PATH="${CONDA_PREFIX}"
    -DCOLVARS_SOURCE_DIR="${colvars_source_dir}"
    -DOPENMM_SOURCE_DIR="${openmm_source_dir}"
    -DLEPTON_DIR="${lepton_dir}"
    -DOPENMM_COLVARS_BUILD_PYTHON=ON
    -DOPENMM_COLVARS_ENABLE_CUDA_KERNELS="${enable_cuda}"
)

if [[ "${enable_cuda}" == "ON" || "${enable_cuda}" == "1" || "${enable_cuda}" == "TRUE" ]]; then
    cmake_args+=(
        -DOPENMM_CUDA_LIBRARY="${openmm_cuda_library}"
        -DOPENMM_COLVARS_CUDA_TOOLKIT_INCLUDE_DIR="${cuda_include_dir}"
    )
fi

echo "Configuring OpenMM-Colvars"
echo "  build dir: ${build_dir}"
echo "  OpenMM env: ${CONDA_PREFIX}"
echo "  OpenMM source: ${openmm_source_dir}"
echo "  Colvars source: ${colvars_source_dir}"
echo "  CUDA kernels: ${enable_cuda}"
if [[ "${enable_cuda}" == "ON" || "${enable_cuda}" == "1" || "${enable_cuda}" == "TRUE" ]]; then
    echo "  OpenMM CUDA library: ${openmm_cuda_library}"
    echo "  CUDA include dir: ${cuda_include_dir}"
fi

cmake "${cmake_args[@]}"

cmake --build "${build_dir}" -j
cmake --install "${build_dir}"

python - <<'PY'
from openmm_colvars import AlchemicalPMEForce, ColvarsForce, ColvarsTABFIntegrator
from openmm_fep import FEPAtoms, FEPConfig, prepare_lambda_abf_system
print("Installed openmm_colvars:", ColvarsForce, AlchemicalPMEForce, ColvarsTABFIntegrator)
print("Installed openmm_fep:", FEPAtoms, FEPConfig, prepare_lambda_abf_system)
PY
