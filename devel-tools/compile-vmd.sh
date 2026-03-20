#!/bin/bash

set -e

source $(dirname $0)/load-recent-git.sh

source $(dirname $0)/set-ccache.sh

# Set explicit path for non-CMake targets
if [ -n "${CCACHE_HOME}" ] ; then
    export PATH=${CCACHE_HOME}:${PATH}
fi

detect_os() {
    # TODO Install dependencies for non RedHat-style OSs?
    if ! declare -f need_rpm 2> /dev/null ; then
        if [ -e /etc/redhat-release ] ; then
            need_rpm() {
                for arg in "${@}" ; do
                    if ! rpm -q -- "${arg}" ; then
                        echo "Error: package ${arg} is not installed." >&2
                        exit 1
                    fi
                done
            }
        else
            need_rpm() { :; }
        fi
    fi
}


check_vmd_dependencies() {

    # Plugins' dependencies
    need_rpm texlive-latex texlive-metafont texlive-pdftex texlive-collection-fontsrecommended

    # VMD dependencies
    need_rpm tcl-devel tk-devel fltk-devel
    need_rpm netcdf-devel
    need_rpm mesa-libGL-devel

    export TCL_VERSION=8.5
    if [ -f /usr/lib64/libtcl8.6.so ] || \
        [ -f /usr/lib/x86_64-linux-gnu/libtcl8.6.so ] ; then
        export TCL_VERSION=8.6
    fi
}


fix_vmd_configure() {

    local configure=$1

    if [ -z "${configure}" ] || [ ! -f "${configure}" ] ; then
        echo "Error: need a path to a VMD configure script" >& 2
        return 1
    fi

    if [ -z "${TCL_VERSION}" ] ; then
        echo "Error: TCL_VERSION could not be defined automatically" >& 2
        return 1
    fi

    sed -i "s/-ltcl8.5/-ltcl${TCL_VERSION}/" ${configure}
    sed -i "s/-ltk8.5/-ltk${TCL_VERSION}/" ${configure}

    sed -i 's/CXX11/CXX17/' ${configure}
    sed -i 's/config_cxx11/config_cxx17/' ${configure}
    sed -i 's/$config_cxx17           = 0/$config_cxx17           = 1/' ${configure}
    sed -i 's/ -fno-for-scope//' ${configure}

    # Also use this to append the flag to nvcc
    sed -i 's/$arch_opt_flag .= " -std=c++11"/$arch_opt_flag .= " -std=c++17"; $arch_nvccflags .= " -std=c++17"/' ${configure}

    # We shouldn't have to install tcsh just for this
    sed -i -e 's/if \[ ! -x "\/bin\/csh/if \[ ! -x "\/bin\/csh-skip-it/' ${configure}

    # VMD2 fixes
    sed -i 's/lfltk_jpeg/ljpeg/' ${configure}
    sed -i 's/lfltk_png/lpng/' ${configure}
    sed -i 's/lfltk_z/lz/' ${configure}

    # Remove ancient CUDA architectures
    sed -i 's/"-gencode arch=compute_30,code=compute_30 "/""/' ${configure}
    sed -i 's/"-gencode arch=compute_30,code=sm_35 "/""/' ${configure}
    sed -i 's/"-gencode arch=compute_30,code=sm_37 "/""/' ${configure}
}


compile_vmd_target() {

    check_vmd_dependencies

    local VMD2PROTOTYPE="no"

    local VMDSRCDIR="${PWD}"
    if [ -f "${1}/src/VMDApp.h" ] ; then
        VMDSRCDIR=$(realpath "${1}")
        shift
    fi
    if [ -f "${1}/vmd/src/VMDApp.h" ] ; then
        VMD2PROTOTYPE="yes"
        VMDSRCDIR=$(realpath "${1}/vmd")
        shift
    fi
    pushd "${VMDSRCDIR}"

    fix_vmd_configure configure

    local label="${1:-multicore}"
    local ret_code

    local VMD_VERSION=$(basename $PWD)

    export VMDINSTALLBINDIR=${1}
    if [ -z "${VMDINSTALLBINDIR}" ] ; then
        if [ "${VMD2PROTOTYPE}" == "yes" ] ; then
            export VMDINSTALLBINDIR="${VMDSRCDIR%/vmd}/install"
        else
            export VMDINSTALLBINDIR="${VMDSRCDIR}/install"
        fi
    fi

    mkdir -p "${VMDINSTALLBINDIR}"

    local VMDINSTALLLIBRARYDIR=${VMDINSTALLBINDIR}

    local -a VMD_OPTS=(TCL TK FLTK OPENGL PTHREADS COLVARS NETCDF)

    VMD_OPTS+=(CXX17)

    export TCL_LIBRARY_DIR=${TCL_LIBRARY_DIR:-/usr/lib64}
    export TCL_INCLUDE_DIR=${TCL_INCLUDE_DIR:-/usr/include}

    export PYTHON_VERSION=$(python3 --version | cut -d' ' -f 2)
    export PYTHON_NAME=python${PYTHON_VERSION%.*}

    if hash python3-config ; then

        export PYTHON_NAME="${PYTHON_NAME}$(python3-config --abiflags)"

        export PYTHON_INCLUDE_DIR="$(python3-config --includes | cut -d' ' -f1)"
        export PYTHON_INCLUDE_DIR=${PYTHON_INCLUDE_DIR#-I}

        export PYTHON_LIB="$(python3-config --libs | cut -d' ' -f1)"
        if [ -z "${PYTHON_LIB}" ] || { echo "${PYTHON_LIB}" | grep -vq lpython; } ; then
            # Not detected, try generating from the Python executable name
            export PYTHON_LIB=-l${PYTHON_NAME}
        fi

        export PYTHON_LIBRARY_DIR="$(python3-config --ldflags | cut -d' ' -f1)"
        export PYTHON_LIBRARY_DIR=${PYTHON_LIBRARY_DIR#-L}
        if [ -z "${PYTHON_LIBRARY_DIR}" ] && [ -d "/usr/lib64" ] ; then
            # Fix for when it's reported empty
            export PYTHON_LIBRARY_DIR="/usr/lib64"
        fi
    fi

    if [ -n "${PYTHON_LIBRARY_DIR}" ] && [ -n "${PYTHON_INCLUDE_DIR}" ] && [ -f "${PYTHON_INCLUDE_DIR}/Python.h" ] ; then

        VMD_OPTS+=(PYTHON)
        sed -i "s/-lpython2.5/${PYTHON_LIB}/" configure

        export PYTHON_PACKAGES_DIR=${PYTHON_LIBRARY_DIR}/${PYTHON_NAME}/site-packages
        if [ ! -d "${PYTHON_PACKAGES_DIR}" ] ; then
            export PYTHON_PACKAGES_DIR=${PYTHON_LIBRARY_DIR}/${PYTHON_NAME%m}/site-packages
        fi

        # NumPy variables
        if [ -n "${VIRTUAL_ENV}" ] ; then
            # Get NumPy from the virtualenv if defined
            export NUMPY_LIBRARY_DIR=${VIRTUAL_ENV}/lib/${PYTHON_NAME}/site-packages/numpy/core
        else
            export NUMPY_LIBRARY_DIR=${NUMPY_LIBRARY_DIR:-${PYTHON_PACKAGES_DIR}/numpy/core}
        fi
        export NUMPY_INCLUDE_DIR=${NUMPY_LIBRARY_DIR}/include

        if [ -n "${NUMPY_LIBRARY_DIR}" ] && [ -d "${NUMPY_LIBRARY_DIR}" ] ; then
            VMD_OPTS+=(NUMPY)
        fi
    fi

    export VMDPLUGINDIR=${VMDINSTALLBINDIR}/plugins

    local VMDPLUGINSRCDIR
    if [ -f ../plugins/Make-arch ] ; then
        # vmd2prototype-like location
        VMDPLUGINSRCDIR=../plugins
    elif [ -f ../plugins/${VMD_VERSION}/Make-arch ] ; then
        VMDPLUGINSRCDIR=../plugins/${VMD_VERSION}
    elif [ -f ../vmd-plugins-source/Make-arch ] ; then
        # Location used in CI tests
        VMDPLUGINSRCDIR=../vmd-plugins-source
    fi

    # Compile and install plugins
    if pushd "${VMDPLUGINSRCDIR}" ; then
        sed -i "s/-ltcl8.5/-ltcl${TCL_VERSION}/" Make-arch
        SQLITEDYNAMIC=1 SQLITEINC=-I/usr/include SQLITELIB=-L/usr/lib64 SQLITELDFLAGS=-lsqlite3 \
            TCLINC=-I/usr/include TCLLIB=-L/usr/lib64 \
            NETCDFINC=-I/usr/include NETCDFLIB=-L/usr/lib64 \
            nice -n 10 make -j$(nproc --all) LINUXAMD64
        PLUGINDIR=${VMDPLUGINDIR} make distrib
        if [ $? = 0 ] ; then
            make clean
        fi
        popd
    else
        echo "Cannot find plugins source repository (tried: ${VMDPLUGINSRCDIR})." >&2
        return 1
    fi

    if [ -n "${CUDA_HOME}" ] && [ -d ${CUDA_HOME} ] ; then
        sed -i 's/\/usr\/local\/cuda-10.2/$ENV{CUDA_HOME}/' configure
        VMD_OPTS+=(CUDA)
        if [ -n "${OPTIX_HOME}" ] && [ -d ${OPTIX_HOME} ] ; then
            sed -i 's/"\/usr\/local\/encap\/NVIDIA-OptiX-SDK-6.5.0-linux64"/$ENV{OPTIX_HOME}/' configure
            VMD_OPTS+=(LIBOPTIX)
        fi
    fi

    if { grep -q LEPTON configure ; } ; then
        VMD_OPTS+=(LEPTON)
    fi

    rm -f plugins
    ln -s ${VMDPLUGINDIR} plugins

    if VMDINSTALLBINDIR=${VMDINSTALLBINDIR} \
        VMDINSTALLLIBRARYDIR=${VMDINSTALLLIBRARYDIR} \
        ./configure LINUXAMD64 ${VMD_OPTS[@]} ; then
        pushd src/
        nice -n 10 make -j$(nproc --all) && make install
        retcode=$?
        if [ $retcode = 0 ] ; then
            make clean
        fi
        popd
    fi

    popd
    return ${retcode}
}


detect_os && compile_vmd_target "${@}"
