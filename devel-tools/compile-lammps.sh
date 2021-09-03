#!/bin/bash


source $(dirname $0)/set-ccache.sh


compile_lammps_target() {

    local CMAKE=cmake
    if hash cmake3 >& /dev/null ; then
        CMAKE=cmake3
    fi
    if hash ${CMAKE} >& /dev/null ; then
        CMAKE_VERSION=$(${CMAKE} --version | head -1 | cut -d' ' -f3)
    else
        echo "Error: no CMake found." >& 2
    fi

    local LAMMPS_SRC_DIR=""
    if [ -f "${1}/src/lmptype.h" ] ; then
        LAMMPS_SRC_DIR=$(realpath "${1}")
        shift
    else
        LAMMPS_SRC_DIR="${PWD}"
    fi

    pushd "${LAMMPS_SRC_DIR}"

    local -a LAMMPS_BUILD_OPTS=()
    local LAMMPS_BUILD_TYPE=Release
    local LAMMPS_INSTALL_DIR="${LAMMPS_SRC_DIR}/install"

    while [ $# -ge 1 ]; do
        if [ "${1,,}" = "kokkos" ] || [ "${1,,}" = "kk" ]; then
            # TODO make the arch selection more portable
            LAMMPS_BUILD_OPTS+=("-DPKG_KOKKOS=on" "-DKokkos_ARCH_HSW=yes" "-DKokkos_ENABLE_OPENMP=yes")
        elif [ "${1,,}" = "debug" ]; then
            LAMMPS_BUILD_TYPE=Debug
            LAMMPS_BUILD_OPTS+=(-DCMAKE_BUILD_TYPE=Debug -DCMAKE_VERBOSE_MAKEFILE=yes)
            LAMMPS_BUILD_OPTS+=(-DCOLVARS_DEBUG)
        else
            LAMMPS_INSTALL_DIR=${1}
        fi
        shift
    done

    if hash mpicc >& /dev/null ; then
        LAMMPS_BUILD_OPTS+=("-DBUILD_MPI=ON")
    fi

    LAMMPS_BUILD_OPTS+=("-DBUILD_OMP=yes")

    # Select FFT
    if [ -n "${CXX}" ] ; then
        if [ $(basename "${CXX}") = "icpc" ] ; then
            LAMMPS_BUILD_OPTS+=("-DFFT=MKL" "-DMKL_INCLUDE_DIRS=${MKLROOT}/include" "-DMKL_LIBRARIES=${MKLROOT}/lib/intel64")
        else
            LAMMPS_BUILD_OPTS+=("-DFFT=FFTW3")
        fi
    else
        LAMMPS_BUILD_OPTS+=("-DFFT=FFTW3")
    fi

    LAMMPS_BUILD_OPTS+=("-DPKG_PYTHON=on")

    if [ -z "${LAMMPS_BUILD_DIR}" ] ; then
        LAMMPS_BUILD_DIR=$(mktemp -d /tmp/lammps-build-XXXXXX)
    fi

    mkdir -p "${LAMMPS_BUILD_DIR}"

    local ret_code=0

    ${CMAKE} \
        -DCMAKE_INSTALL_PREFIX="${LAMMPS_INSTALL_DIR}" \
        -DCMAKE_BUILD_TYPE=${LAMMPS_BUILD_TYPE:-Release} \
        -DBUILD_SHARED_LIBS=on \
        -C "${LAMMPS_SRC_DIR}/cmake/presets/most.cmake" \
        -C "${LAMMPS_SRC_DIR}/cmake/presets/nolib.cmake" \
        ${LAMMPS_BUILD_OPTS[@]} \
        -S "${LAMMPS_SRC_DIR}/cmake" \
        -B "${LAMMPS_BUILD_DIR}" \
        && \
        ${CMAKE} --build "${LAMMPS_BUILD_DIR}" --parallel $(nproc --all)
    ret_code=$?

    if [ $ret_code = 0 ] && [ -n "${LAMMPS_INSTALL_DIR}" ] ; then
        pushd "${LAMMPS_BUILD_DIR}"
        make install
        ret_code=$?
        if [ $ret_code = 0 ] ; then
            rm -fr "${LAMMPS_BUILD_DIR}"
        fi
        popd
    fi

    popd

    return ${ret_code}
}


compile_lammps_target "${@}"
