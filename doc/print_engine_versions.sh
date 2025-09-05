#!/bin/bash

TOPDIR=$(git rev-parse --show-toplevel)
if [ ! -d ${TOPDIR} ] ; then
  echo "Error: cannot identify top project directory." >& 2
  exit 1
fi

source ${TOPDIR}/devel-tools/version_functions.sh

if [ -n "${GITHUB_ACTION}" ] ; then
    export remote=origin
fi

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

    local -a branches
    if [ -n "${remote}" ] ; then
        branches=($(git branch -r --format='%(refname)' | grep "refs/remotes/${remote}" | sed -s "s/refs\/remotes\/${remote}\///" | grep ^${tag_prefix}))
    else
        branches=($(git branch -l --format='%(refname)' | sed -s 's/refs\/heads\///' | grep ^${tag_prefix}))
    fi

    for branch in ${branches[@]} ; do
        local package_version=${branch#${tag_prefix}}
        if [ -n "${remote}" ] ; then
            colvars_version=$(get_colvarmodule_version ${remote}/${branch})
        else
            colvars_version=$(get_colvarmodule_version ${branch})
        fi
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

echo "This page only lists major (feature) releases of simulation packages, not patch releases (e.g. NAMD 3.0.2). Patch release generally do not contain a specific version of the Colvars master branch, but the same version as the previous major release, possibly with later fixes added."

for package in GROMACS LAMMPS NAMD VMD ; do
    echo -n " [[${package}](#versions-included-in-${package})]"
done
echo
echo


for package in GROMACS LAMMPS NAMD VMD ; do
    echo "### Versions included in ${package}"
    echo "${package} version | Colvars version"
    echo "-------------- | ---------------"
    sort_command=sort_versions
    if [ ${package} == LAMMPS ] ; then
        sort_command=sort_lammps_versions
    fi
    print_tag_versions ${package} ${package,,}- | ${sort_command}
    echo
done

echo
echo "This page was last updated on: $(date +'%Y-%m-%d')"
