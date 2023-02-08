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

# Create a directory (or find an existing one) and set the variable $dirname
create_test_dir() {

    local workdir=${WORKDIR%/}

    if ls $workdir/ | grep -q "^[0-9]\{3\}_${1}$" ; then
        dirname=`ls $WORKDIR/ | grep "^[0-9]\{3\}_${1}$"`
        dirname=${WORKDIR}/${dirname}
        return
    fi

    dirname="${workdir}/"`printf %03d ${n_test}`_$1
    mkdir -p ${dirname}/AutoDiff && \
        echo "Created new directory $dirname"
}


write_colvars_config() {
    local colvar=$1
    local bias=$2
    local filename=$3
    if [ -z "${filename}" ] ; then
        filename=${dirname}/test.in
    fi
    echo 'colvarsTrajFrequency 1' > ${filename}
    echo 'colvarsRestartFrequency 10' >> ${filename}
    cat indexfile.in >> ${filename}
    echo '' >> ${filename}
    cat ${colvar}.in >> ${filename}
    if [ "x$bias" != "x" ] ; then
        echo '' >> ${filename}
        cat ${bias}.in >> ${filename}
    fi
}

dirname=''
n_test=0

m4 < harmonic.in.m4 > harmonic-fixed.in
m4 -Dcenters_moving_half < harmonic.in.m4 > harmonic-centers-moving.in
m4 -Dcenters_moving_full < harmonic.in.m4 > harmonic-centers-moving-full.in
m4 -Dcenters_moving_stages < harmonic.in.m4 > harmonic-centers-moving-stages.in


m4 -Dti_pmf < harmonic.in.m4 > harmonic-fixed-ti.in
m4 -Dti_pmf -Dcenters_moving_half < harmonic.in.m4 > harmonic-centers-moving-ti.in


for colvar in "distance" ; do
    for bias in \
        "harmonic-fixed" \
        "harmonic-centers-moving" \
        "harmonic-centers-moving-full" \
        "harmonic-centers-moving-stages" \
        "harmonic-k-moving" \
        "harmonic-k-moving-stages" \
        "harmonicwalls-upper-fixed" \
        "harmonicwalls-both-fixed" \
        "harmonicwalls-k-moving" \
        "harmonicwalls-upper-k-moving" \
        "harmonicwalls-upper-k-moving-decoupling" \
        "harmonicwalls-both-k-moving" \
        "linear-fixed" \
        "linear-k-moving" \
        ; do
        create_test_dir ${colvar}_${bias}
        write_colvars_config ${colvar} ${bias}
    done
done

create_test_dir "distance_histogram-customgrid"
write_colvars_config "distance" "histogram-customgrid"

create_test_dir "dihedral_histogram-customgrid-periodic"
write_colvars_config "dihedral" "histogram-customgrid-periodic"

create_test_dir "distance-dihedral_histogram2d-periodic"
write_colvars_config "distance-dihedral_histogram2d-periodic"

create_test_dir "distance-dihedral_histogram2d-periodic-customgrid"
write_colvars_config "distance-dihedral_histogram2d-periodic-customgrid"

create_test_dir "distancewalls"
write_colvars_config "distance" "harmonicwalls-fixed"
write_colvars_config "distancewalls" "" ${dirname}/test.legacy.in

create_test_dir "distancewalls-compatible"
write_colvars_config "distancewalls" ""

create_test_dir "dihedralwalls"
write_colvars_config "dihedral" "harmonicwalls360angle-fixed"
write_colvars_config "dihedralwalls" "" ${dirname}/test.legacy.in

create_test_dir "distance-wall-bypassExtended-off"
write_colvars_config "distance-extended" "harmonicwalls-bypassExtended-off"

# Tests for each colvar without a bias
for colvar in \
    "distance" \
    ; do
    create_test_dir ${colvar}
    write_colvars_config ${colvar} ""
    # if [ -f ${colvar}-fitgroup.in ] ; then
    #     create_test_dir ${colvar}-fitgroup_${bias}
    #     write_colvars_config ${colvar}-fitgroup ${bias}
    #     m4 -DfittingGroup=refPositionsGroup < ${colvar}-fitgroup.in > ${colvar}-refposgroup.in
    #     write_colvars_config ${colvar}-refposgroup ${bias} ${dirname}/test.legacy.in
    # fi
done

create_test_dir "distance-extended"
write_colvars_config "distance-extended" ""

create_test_dir "distance-runave"
write_colvars_config "distance-runave" ""

create_test_dir "distance-autocorrfunc"
write_colvars_config "distance-autocorrfunc" ""

create_test_dir "distance-corrfunc"
write_colvars_config "distance-corrfunc" ""


for colvar in "distance-grid" ; do
    for bias in \
        "harmonic-fixed" \
        "harmonic-centers-moving" \
        "harmonic-fixed-ti" \
        "harmonic-centers-moving-ti" \
        "harmonic-k-moving" \
        "histogram" \
        "histogram-customgrid" \
        "abf" \
        "metadynamics" \
        "metadynamics-sigmas" \
        "metadynamics-wt" \
        "metadynamics-ti" \
        ; do
        create_test_dir ${colvar}_${bias}
        write_colvars_config ${colvar} ${bias}
    done
done

create_test_dir "distance-grid-expand_metadynamics"
write_colvars_config "distance-grid-expand" "metadynamics"

create_test_dir "distance-grid_metadynamics-rebin"
write_colvars_config "distance-grid" "metadynamics-keephills"
write_colvars_config "distance-grid-finer" "metadynamics-keephills-rebingrids" ${dirname}/test.restart.in

create_test_dir "distance-grid_metadynamics-rebinfromgrids"
write_colvars_config "distance-grid" "metadynamics-new"
write_colvars_config "distance-grid-finer" "metadynamics-rebingrids" ${dirname}/test.restart.in

colvar="distance-extended-grid"
for bias in "histogram" "histogram-bypassExtended"; do
    create_test_dir ${colvar}_${bias}
    write_colvars_config ${colvar} ${bias}
done

m4 -Daxis < distancez.in.m4 > distancez-axis.in
m4 -Daxis -Dfitgroup < distancez.in.m4 > distancez-axis-fitgroup.in
m4 < distancez.in.m4 > distancez.in
m4 -Dfitgroup < distancez.in.m4 > distancez-fitgroup.in
m4 -Daxis -DdistanceZ=distanceXY < distancez.in.m4 > distancexy-axis.in

m4 -Dorientation=orientationAngle < orientation.in > orientationangle.in
m4 -Dorientation=orientationAngle < orientation-fitgroup.in > orientationangle-fitgroup.in
m4 -Dorientation=orientationProj < orientation.in > orientationproj.in
m4 -Dorientation=orientationProj < orientation-fitgroup.in > orientationproj-fitgroup.in

# Tests of individual collective variables
for colvar in \
    "angle" \
    "dihedral" \
    "coordnum" "coordnum-group2centeronly" "coordnum-pairlist" \
    "coordnum-aniso" "coordnum-aniso-pairlist" \
    "groupcoord" "groupcoord-aniso" \
    "dipoleangle" "dipolemagnitude" \
    "distancez" "distancez-fitgroup" \
    "distancez-axis" "distancez-axis-fitgroup" \
    "distancexy-axis" \
    "distanceinv" \
    "distance-coeffs" \
    "gyration" \
    "hbond" \
    "inertia" \
    "inertiaz" \
    "rmsd" \
    "eigenvector" \
    "eigenvector-difference" \
    "eigenvector-normalized" \
    "tilt" \
    "spinangle" \
    "selfcoordnum" "selfcoordnum-pairlist" \
    "orientationangle" "orientationproj" \
   ; do
    for bias in \
        "harmonic-fixed" \
        ; do
        create_test_dir ${colvar}_${bias}
        write_colvars_config ${colvar} ${bias}
        if [ -f ${colvar}-fitgroup.in ] ; then
            create_test_dir ${colvar}-fitgroup_${bias}
            write_colvars_config ${colvar}-fitgroup ${bias}
            m4 -DfittingGroup=refPositionsGroup < ${colvar}-fitgroup.in > ${colvar}-refposgroup.in
            write_colvars_config ${colvar}-refposgroup ${bias} ${dirname}/test.legacy.in
        fi
    done
done

colvar="orientation"
for bias in "harmonic-ori-fixed" "harmonic-ori-moving" ; do
create_test_dir ${colvar}_${bias}
write_colvars_config ${colvar} ${bias}
if [ -f ${colvar}-fitgroup.in ] ; then
    create_test_dir ${colvar}-fitgroup_${bias}
    write_colvars_config ${colvar}-fitgroup ${bias}
    m4 -DfittingGroup=refPositionsGroup < ${colvar}-fitgroup.in > ${colvar}-refposgroup.in
    write_colvars_config ${colvar}-refposgroup ${bias} ${dirname}/test.legacy.in
fi
done


create_test_dir customfunction_harmonic-fixed
write_colvars_config customfunction_harmonic-fixed ""


colvar="distancevec"
bias="harmonic-dvec-fixed"
create_test_dir ${colvar}_${bias}
write_colvars_config ${colvar} ${bias}
if [ -f ${colvar}-fitgroup.in ] ; then
    create_test_dir ${colvar}-fitgroup_${bias}
    write_colvars_config ${colvar}-fitgroup ${bias}
    m4 -DfittingGroup=refPositionsGroup < ${colvar}-fitgroup.in > ${colvar}-refposgroup.in
    write_colvars_config ${colvar}-refposgroup ${bias} ${dirname}/test.legacy.in
fi

colvar="distancedir"
for bias in "harmonic-ddir-fixed" "harmonic-ddir-moving" ; do
    create_test_dir ${colvar}_${bias}
    write_colvars_config ${colvar} ${bias}
done
bias="harmonic-ddir-fixed"
create_test_dir ${colvar}-fitgroup_${bias}
write_colvars_config ${colvar}-fitgroup ${bias}


colvar="dihedralPC"
bias="abf2d"
create_test_dir ${colvar}_${bias}
write_colvars_config ${colvar} ${bias}

colvar="dihedralPC"
bias="metadynamics-2d"
create_test_dir ${colvar}_${bias}
write_colvars_config ${colvar} ${bias}

colvar="distancepairs"
bias="linear-distancepairs"
create_test_dir ${colvar}_${bias}
write_colvars_config ${colvar} ${bias}

colvar="distancepairs"
bias="histogramrestraint-dp"
create_test_dir ${colvar}_${bias}
write_colvars_config ${colvar} ${bias}

colvar="eulerangles"
bias="harmonic-fixed-euler"
create_test_dir ${colvar}_${bias}
write_colvars_config ${colvar} ${bias}


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

