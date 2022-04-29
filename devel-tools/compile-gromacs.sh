#!/bin/bash

source $(dirname $0)/load-recent-git.sh

source $(dirname $0)/load-recent-gcc.sh

source $(dirname $0)/set-ccache.sh


compile_gromacs_target() {

    local CMAKE=cmake CTEST=ctest
    if hash cmake3 >& /dev/null ; then
        CMAKE=cmake3
        CTEST=ctest3
    fi
    if hash ${CMAKE} >& /dev/null ; then
        CMAKE_VERSION=$(${CMAKE} --version | head -1 | cut -d' ' -f3)
    else
        echo "Error: no CMake found." >& 2
    fi

    local GMX_SRC_DIR=""
    if [ -f "${1}/src/gromacs/commandline.h" ] ; then
        GMX_SRC_DIR=$(realpath "${1}")
        shift
    else
        GMX_SRC_DIR="${PWD}"
    fi

    pushd "${GMX_SRC_DIR}"

    local -a GMX_BUILD_OPTS=()
    local GMX_BUILD_TYPE=Release
    local GMX_INSTALL_DIR="${GMX_SRC_DIR}/install"

    if [ -n "${GITHUB_ACTION}" ] ; then
        # We don't want to let the "dirty" flag constantly invalidate the cache
        GMX_BUILD_OPTS+=("-DGMX_GIT_VERSION_INFO=OFF")
    fi

    GMX_BUILD_OPTS+=("-DGMX_DOUBLE=ON")

    while [ $# -ge 1 ]; do
        if [ "${1,,}" = "debug" ]; then
            GMX_BUILD_TYPE=Debug
            GMX_BUILD_OPTS+=(-DCMAKE_BUILD_TYPE=Debug -DCMAKE_VERBOSE_MAKEFILE=yes)
            GMX_BUILD_OPTS+=(-DCOLVARS_DEBUG)
        else
            GMX_INSTALL_DIR=${1}
        fi
        shift
    done

    if hash mpicc >& /dev/null ; then
        GMX_BUILD_OPTS+=("-DGMX_MPI=ON")
    fi

    # GMX_BUILD_OPTS+=("-DBUILD_OMP=yes")

    # Select FFT
    if [ "$(basename \"${CXX}\")" = "icpc" ] ; then
        GMX_BUILD_OPTS+=("-DGMX_FFT_LIBRARY=mkl")
    else
        GMX_BUILD_OPTS+=("-DGMX_FFT_LIBRARY=fftw3")
        if [ ! -f /usr/include/fftw3.h ] ; then
            GMX_BUILD_OPTS+=("-DGMX_BUILD_OWN_FFTW=ON")
        fi
    fi

    # When on GitHub Actions, download the tests as well
    if [ -n "${GITHUB_ACTION}" ] ; then
        GMX_BUILD_OPTS+=(-DREGRESSIONTEST_DOWNLOAD=ON)
    fi

    if [ -z "${GMX_BUILD_DIR}" ] ; then
        GMX_BUILD_DIR="${GMX_SRC_DIR}/build"
    fi

    mkdir -p "${GMX_BUILD_DIR}"

    local ret_code=0

    ${CMAKE} \
        -DCMAKE_INSTALL_PREFIX="${GMX_INSTALL_DIR}" \
        -DCMAKE_BUILD_TYPE=${GMX_BUILD_TYPE:-Release} \
        ${GMX_BUILD_OPTS[@]} \
        -S "${GMX_SRC_DIR}" \
        -B "${GMX_BUILD_DIR}" \
        && \
        ${CMAKE} --build "${GMX_BUILD_DIR}" --parallel $(nproc --all)
    ret_code=$?

    if [ -n "${GITHUB_ACTION}" ] ; then
        # On GitHub Actions, build the tests as well
        ${CMAKE} --build "${GMX_BUILD_DIR}" --target tests --parallel $(nproc --all)
        ret_code=$((${ret_code} || $?))
    fi

    if [ ${ret_code} = 0 ] ; then
        pushd "${GMX_BUILD_DIR}"
        ${CTEST} --output-on-failure
        retcode=$?
        if [ -n "${GMX_INSTALL_DIR}" ] ; then
            make install
            ret_code=$((${ret_code} || $?))
        fi
        popd
    fi

    popd

    return ${ret_code}
}


compile_gromacs_target "${@}"
