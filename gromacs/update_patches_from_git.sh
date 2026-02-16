#!/bin/bash

GROMACS_SOURCE=${1}
if [ -n "${GROMACS_SOURCE}" ] ; then
    GROMACS_SOURCE=$(git -C ${GROMACS_SOURCE} rev-parse --show-toplevel)
fi

if [ $? != 0 ] || [ ! -s ${GROMACS_SOURCE}/src/gromacs/version.h.cmakein ] ; then
    echo >&2
    echo "Usage: $0 <path_to_GROMACS_tree>" >& 2
    exit 1
fi

COLVARS_SOURCE=$(dirname $(realpath $0))
COLVARS_SOURCE=$(git -C ${COLVARS_SOURCE} rev-parse --show-toplevel)

modified_files=($(git -C ${GROMACS_SOURCE} ls-files -m src/gromacs/applied_forces/colvars))

for file in ${modified_files[@]} ; do
    git -C ${GROMACS_SOURCE} diff ${file} > ${COLVARS_SOURCE}/gromacs/src/${file#src\/gromacs\/}.patch
done
