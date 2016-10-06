#!/bin/bash

# Create configuration files for Colvars regression tests that only use features shared between codes.
# Colvars config files are copied by this script in a code-specific folder working directory (e.g. namd/tests/library).

WORKDIR=${PWD}
if [ -n "${1}" ] ; then
    if [ -d ${1} ] ; then
        WORKDIR=${1}
    else
        echo "Directory ${1} does not exist"
        exit 1
    fi
fi

create_test_dir() {

    if ls $WORKDIR/ | grep -q "_${1}" ; then
        dirname=`ls $WORKDIR/ | grep "_${1}"`
        echo "$dirname exists"
        dirname=${WORKDIR}/${dirname}
        return
    fi
    
    # Add to an existing set of regression tests
    while ls $WORKDIR/ | grep -q `printf %03d $n_test`_ ; do
        n_test=$((++n_test))
    done

    dirname="${WORKDIR}/`printf %03d ${n_test}`_$1"
    if [ ! -d ${dirname} ] ; then
        mkdir ${dirname}
    fi
}

dirname=''
n_test=1

for colvar in "distance" ; do
    for bias in \
        "harmonic-fixed" \
            "harmonic-centers-moving" \
            "harmonic-k-moving" \
        ; do
        create_test_dir ${colvar}_${bias}
        echo 'colvarsTrajFrequency 1' > ${dirname}/test.in 
        echo 'colvarsRestartFrequency 10' >> ${dirname}/test.in 
        cat indexfile.in >> ${dirname}/test.in
        echo '' >> ${dirname}/test.in 
        cat ${colvar}.in >> ${dirname}/test.in
        echo '' >> ${dirname}/test.in 
        cat ${bias}.in >> ${dirname}/test.in
    done
done

for colvar in "distance-grid" ; do
    for bias in \
        "harmonic-fixed" \
            "harmonic-centers-moving" \
            "harmonic-k-moving" \
            "abf" \
            "metadynamics" \
        ; do
        create_test_dir ${colvar}_${bias}
        echo 'colvarsTrajFrequency 1' > ${dirname}/test.in 
        echo 'colvarsRestartFrequency 10' >> ${dirname}/test.in 
        cat indexfile.in >> ${dirname}/test.in
        echo '' >> ${dirname}/test.in 
        cat ${colvar}.in >> ${dirname}/test.in
        echo '' >> ${dirname}/test.in 
        cat ${bias}.in >> ${dirname}/test.in
    done
done

# TODO uncomment this and the add two-dimensional regtests
# # Generate two-variables versions of bias configurations
# for bias in "harmonic-fixed" \
#     "harmonic-centers-moving" \
#     "harmonic-k-moving" \
#     "abf" \
#     "metadynamics" \
#     ; do
#     cat ${bias}.in \
#         | sed 's/colvars        one/colvars        one two/' \
#         | sed 's/centers        0.0/centers        0.0 0.0/' \
#         | sed 's/targetCenters  0.1/targetCenters  0.1 0.1/' \
#         > ${bias}-2.in
# done

unset -f n_test dirname
