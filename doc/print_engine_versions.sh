#!/bin/bash

source ../devel-tools/version_functions.sh


reformat_lammps_version() {
    local lammpsversion=$1
    python3 -c "
import datetime
print(datetime.datetime.strptime('${lammpsversion}','%d%b%Y').strftime('%Y-%m-%d'))
"
}


print_tag_versions() {
    local package=$1
    local tag_prefix=$2
    local tag
    local colvars_version
    local lammps_version=""
    for tag in $(git tag -l|grep ^${tag_prefix}); do
        local package_version=${tag#${tag_prefix}}
        colvars_version=$(get_colvarmodule_version ${tag})
        if [ ${package} = LAMMPS ] ; then
            lammps_version=$(reformat_lammps_version ${package_version#*_})" "
        fi
        echo "${lammps_version}${package_version} | [${colvars_version}](https://github.com/Colvars/colvars/releases/tag/${tag})"
    done
}

echo "## List of Colvars versions included in simulation/analysis packages"
echo


sort_versions(){
    sort -r
}


sort_lammps_versions(){
    sort -r | cut -d' ' -f2-
}


for package in LAMMPS NAMD VMD ; do
    echo "### Versions included in ${package}"
    echo "${package} version | Colvars version"
    echo "-------------- | ---------------"
    sort_command=sort_versions
    if [ ${package} = LAMMPS ] ; then
        sort_command=sort_lammps_versions
    fi
    print_tag_versions ${package} $(echo ${package}- | tr '[:upper:]' '[:lower:]') | ${sort_command}
    echo
done
