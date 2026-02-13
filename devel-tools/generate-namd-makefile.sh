#!/bin/bash

set -e

# Generate a file suitable for inclusion by the NAMD Makefile

if [ -f "${2}" ] ; then
    echo "Usage: ${0} <sources_var_name> <sources_path> <objs_var_name> <obj_path> [... source files ...]" >& 2
    exit 1
fi

SRC_VAR_NAME=${1}
shift
SRC_PREFIX=${1}
shift
OBJ_VAR_NAME=${1}
shift
OBJ_PREFIX=${1}
shift

files=($(echo "${@}" | sort))

echo "${SRC_VAR_NAME} = \\"
for file in "${files[@]}" ; do
    if [ ! -f "${file}" ] ; then
        echo "Error: file ${file} is missing" >&2
        echo "Usage: ${0} <variable_name> <destination_path> [... source files ...]" >&2
        exit 1
    fi
    file=$(basename ${file})
    echo "	${SRC_PREFIX}/${file} \\"
done
echo ""
echo ""


echo "${OBJ_VAR_NAME} = \\"
for file in "${files[@]}" ; do
    file=$(basename ${file})
    prefix=${file%.C}
    prefix=${prefix%.cpp}
    prefix=${prefix%.cu}
    echo "	${OBJ_PREFIX}/${prefix}.o \\"
done
echo ""
echo ""
