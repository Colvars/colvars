#!/bin/bash

# TODO most of the library paths are for RH-style Linux environments

source $(dirname $0)/set-ccache.sh

detect_os() {

    # if ! hash lsb_release ; then
    #     echo "Error: lsb_release is not installed." >&2
    #     if [ -n "$BASH_SOURCE" ] ; then
    #         return 1
    #     fi
    # fi

    # DIST_NAME=$(lsb_release -is)
    # DIST_VERSION=$(lsb_release -sr | sed -e 's/\..*//')

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


get_first_word() {
    local line="$*"
    line=$(echo ${line})
    echo $(echo ${line} | cut -d' ' -f 1)
}


compile_vmd_target() {

    local VMDSRCDIR="${PWD}"
    if [ -f "${1}/src/VMDApp.h" ] ; then
        VMDSRCDIR=$(realpath "${1}")
        shift
    fi
    pushd "${VMDSRCDIR}"

    local label="${1:-multicore}"
    local ret_code

    local VMD_VERSION=$(basename $PWD)

    export VMDINSTALLBINDIR=${1}
    if [ -z "${VMDINSTALLBINDIR}" ] ; then
        export VMDINSTALLBINDIR="${VMDSRCDIR}/install"
    fi

    mkdir -p "${VMDINSTALLBINDIR}"

    check_vmd_dependencies

    local VMDINSTALLLIBRARYDIR=${VMDINSTALLBINDIR}

    local -a VMD_OPTS=(TCL TK FLTK OPENGL PTHREADS COLVARS NETCDF)

    export TCL_LIBRARY_DIR=${TCL_LIBRARY_DIR:-/usr/lib64}
    export TCL_INCLUDE_DIR=${TCL_INCLUDE_DIR:-/usr/include}

    sed -i "s/-ltcl8.5/-ltcl${TCL_VERSION}/" configure
    sed -i "s/-ltk8.5/-ltk${TCL_VERSION}/" configure

    export PYTHON_NAME=$(basename $(which python3))

    if hash python3-config ; then
        export PYTHON_LIB="$(python3-config --libs | cut -d' ' -f1)"
        export PYTHON_NAME="${PYTHON_NAME}$(python3-config --abiflags)"

        if [ -z "${PYTHON_INCLUDE_DIR}" ] ; then
            export PYTHON_INCLUDE_DIR="$(python3-config --includes | cut -d' ' -f1)"
            export PYTHON_INCLUDE_DIR=${PYTHON_INCLUDE_DIR#-I}
            export PYTHON_LIBRARY_DIR="$(python3-config --ldflags | cut -d' ' -f1)"
            export PYTHON_LIBRARY_DIR=${PYTHON_LIBRARY_DIR#-L}
        fi
    else
        export PYTHON_INCLUDE_DIR=${PYTHON_INCLUDE_DIR:-/usr/include/${PYTHON_NAME}}
        export PYTHON_LIBRARY_DIR=${PYTHON_LIBRARY_DIR:-/usr/lib64}
    fi

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

    if [ -n "${PYTHON_LIBRARY_DIR}" ] && [ -d "${PYTHON_LIBRARY_DIR}" ] ; then
        VMD_OPTS+=(PYTHON)
        sed -i "s/-lpython2.5/${PYTHON_LIB}/" configure
    fi

    if [ -n "${NUMPY_LIBRARY_DIR}" ] && [ -d "${NUMPY_LIBRARY_DIR}" ] ; then
        VMD_OPTS+=(NUMPY)
    fi

    export VMDPLUGINDIR=${VMDINSTALLBINDIR}/plugins

    local VMDPLUGINSRCDIR
    if [ -d ../plugins/${VMD_VERSION} ] ; then
        VMDPLUGINSRCDIR=../plugins/${VMD_VERSION}
    fi
    if [ -d ../vmd-plugins-source ] ; then
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

    if { ! grep -q csh-skip-it configure ; } ; then
        # We should't have to install tcsh just for this
        sed -i -e 's/if \[ ! -x "\/bin\/csh/if \[ ! -x "\/bin\/csh-skip-it/' configure
    fi

    if [ -n "${CUDA_HOME}" ] && [ -d ${CUDA_HOME} ] ; then
        sed -i -e "s/\/usr\/local\/cuda-10.2/${CUDA_HOME////\/}/" configure
        VMD_OPTS+=(CUDA)
        if [ -n "${CUDA_HOME}" ] && [ -d ${CUDA_HOME} ] ; then
            sed -i -e "s/\/usr\/local\/encap\/NVIDIA-OptiX/${OPTIX_HOME////\/}/" configure
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

    if [ $retcode = 0 ] ; then
        rm -f plugins
        unset TCL_LIBRARY_DIR TCL_INCLUDE_DIR PYTHON_LIB \
            PYTHON_NAME PYTHON_INCLUDE_DIR PYTHON_LIBRARY_DIR PYTHON_PACKAGES_DIR \
            NUMPY_INCLUDE_DIR NUMPY_LIBRARY_DIR
    fi

    popd
    return ${retcode}
}


detect_os

compile_vmd_target "${@}"
