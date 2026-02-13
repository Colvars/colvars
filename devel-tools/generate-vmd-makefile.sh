#!/bin/bash

set -e

# Generate the file colvars_files.pl used by VMD (i.e. a poor person's version of a dependency file)

if [ ! -f "${1}" ] ; then
    echo "Usage: ${0} [... source files ...] [... header files ...]" >& 2
    exit 1
fi

all_files=("${@}")


canonize_filename() {
    local file=$1
    if [ "${file%.cpp}" != "${file}" ] ; then
        basename "${file%.cpp}.C"
    else
        basename "${file}"
    fi
}

source_files=()
header_files=()
cuda_files=()
for file in "${all_files[@]}" ; do
    if [ "${file%.cpp}" != "${file}" ] ; then
        source_files+=($(canonize_filename "${file}"))
    elif [ "${file%.C}" != "${file}" ] ; then
        source_files+=($(canonize_filename "${file}"))
    elif [ "${file%.h}" != "${file}" ] ; then
        header_files+=($(canonize_filename "${file}"))
    elif [ "${file%.cu}" != "${file}" ] ; then
        echo "Error: building of Colvars CUDA files with VMD is not supported yet" >&2
        # cuda_files+=($(canonize_filename "${file}"))
    else
        echo "Error: file ${file} has an unsupported extension" >&2
        exit 1
    fi
done

cat <<EOF
# List of files for the Colvars module

our ( \$colvars_defines );

our ( @colvars_cc );
our ( @colvars_cu );
our ( @colvars_ccpp );
our ( @colvars_h );

\$colvars_defines = " -DVMDCOLVARS";

@colvars_cc = ();
EOF


print_list() {
    local var_name=$1
    shift
    local -a files=($(echo "${@}" | tr ' ' '\n' | sort | xargs))
    echo "@${var_name} = ("
    for file in "${files[@]:0:${#files[@]}-1}" ; do
        echo "    '${file}',"
    done
    echo "    '${files[${#files[@]}-1]}'"
    echo "    );"
}

if [ -n "${cuda_files}" ] ; then
    print_list "colvars_cu" "${cuda_files[@]}"
else
    echo "@colvars_cu = ();"
fi

if [ -n "${source_files}" ] ; then
    print_list "colvars_ccpp" "${source_files[@]}"
else
    echo "@colvars_ccpp = ();"
fi

if [ -n "${header_files}" ] ; then
    print_list "colvars_h" "${header_files[@]}"
else
    echo "@colvars_h = ();"
fi
