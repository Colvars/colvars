#!/bin/bash

TOPDIR=$(git rev-parse --show-toplevel)
if [ ! -d ${TOPDIR} ] ; then
    echo "Error: cannot identify top project directory." >& 2
    exit 1
fi

tests=([0-9][0-9][0-9]_*)
tests=(${tests[@]/010_MPI_multi-sim})

../library/run_tests.sh ${tests[@]}

# Run tests that depend on an MPI build
if source ${TOPDIR}/devel-tools/load-openmpi.sh ; then
    if cd 010_MPI_multi-sim ; then
        ./run.sh
        cd -
    fi
fi
