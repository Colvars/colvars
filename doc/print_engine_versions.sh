#!/bin/bash

TOPDIR=$(git rev-parse --show-toplevel)
if [ ! -d ${TOPDIR} ] ; then
  echo "Error: cannot identify top project directory." >& 2
  exit 1
fi

source ${TOPDIR}/devel-tools/version_functions.sh


reformat_lammps_version() {
    local lammpsversion=$1
    if [ ${lammpsversion} == develop ] ; then
        echo ${lammpsversion}
    else
        python3 -c "
import datetime
print(datetime.datetime.strptime('${lammpsversion}','%d%b%Y').strftime('%Y-%m-%d'))
"
    fi
}


print_tag_versions() {
    local package=$1
    local tag_prefix=$2
    local tag
    local colvars_version
    local lammps_version=""

    for tag in $(git tag -l | grep ^${tag_prefix}) ; do
        local package_version=${tag#${tag_prefix}}
        colvars_version=$(get_colvarmodule_version ${tag})
        if [ ${package} == LAMMPS ] ; then
            lammps_version=$(reformat_lammps_version ${package_version#*_})" "
        fi
        echo "${lammps_version}${package_version} | [${colvars_version}](https://github.com/Colvars/colvars/releases/tag/${tag})"
    done

    for branch in $(git branch -l --format='%(refname)' | sed -s 's/refs\/heads\///' | grep ^${tag_prefix}) ; do
        local package_version=${branch#${tag_prefix}}
        colvars_version=$(get_colvarmodule_version ${branch})
        if [ ${package} == LAMMPS ] ; then
            lammps_version=$(reformat_lammps_version ${package_version#*_})" "
        fi
        if [ ${package} == GROMACS ] ; then
            if [ ${branch} == gromacs-2023 ] || [ ${branch} == gromacs-2022 ] ; then
                # These branches do not reflect GROMACS standard releases
                continue
            fi
        fi
        echo "${lammps_version}${package_version} | [${colvars_version}](https://github.com/Colvars/colvars/tree/${branch})"
    done
}


sort_versions(){
    sort -r
}


sort_lammps_versions(){
    sort -r | cut -d' ' -f2-
}


for package in GROMACS LAMMPS NAMD VMD ; do
    echo -n " [[${package}](#versions-included-in-${package})]"
done
echo
echo


for package in GROMACS LAMMPS NAMD VMD ; do
    echo "### Versions included in ${package}"
    if [ ${package} == NAMD ] ; then
        echo "(Note: the Colvars version included in NAMD 2.12 is the same as the one included in 2.12b1 with only bugfixes applied: therefore, NAMD 2.12 does not correspond to a specific version of the Colvars source tree)"
        echo
    fi
    echo "${package} version | Colvars version"
    echo "-------------- | ---------------"
    sort_command=sort_versions
    if [ ${package} == LAMMPS ] ; then
        sort_command=sort_lammps_versions
    fi
    print_tag_versions ${package} ${package,,}- | ${sort_command}
    echo
done
