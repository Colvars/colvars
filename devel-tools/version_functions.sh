# -*- bash -*-
# Functions to manage version strings


get_colvars_macro() {
    local commit="${1}:"
    if [ "${commit}" = ":" ] ; then
        # Use work tree
        commit=""
    fi
    local macro=$2
    local file1=$3
    local file2=$4
    local version=$(git grep "#define ${macro}" ${commit}${file1} 2> /dev/null)
    if [ "x${version}" = "x" ] ; then
        version=$(git grep "#define ${macro}" ${commit}${file2} 2> /dev/null)
    fi
    version=$(echo ${version} | awk '{ print $3 }' | tr -d '"')
    echo ${version}
}

get_colvarmodule_version() {
    local commit=$1
    get_colvars_macro "${commit}" COLVARS_VERSION \
                      src/colvars_version.h src/colvarmodule.h
}

get_colvarproxy_lammps_version() {
    local commit=$1
    get_colvars_macro "${commit}" COLVARPROXY_VERSION \
                      lammps/src/USER-COLVARS/colvarproxy_lammps_version.h \
                      lammps/src/USER-COLVARS/colvarproxy_lammps.h
}

get_colvarproxy_namd_version() {
    local commit=$1
    get_colvars_macro "${commit}" COLVARPROXY_VERSION \
                      namd/src/colvarproxy_namd_version.h \
                      namd/src/colvarproxy_namd.h
}

get_colvarproxy_vmd_version() {
    local commit=$1
    get_colvars_macro "${commit}" COLVARPROXY_VERSION \
                      vmd/src/colvarproxy_vmd_version.h \
                      vmd/src/colvarproxy_vmd.h
}

get_last_master_commit() {
    local grep_pattern=${1}
    local commit
    local last_master

    local first_since_master=$(git log --pretty=format:'%h' \
                                   --date=short master..HEAD | tail -n 1)
    if [ "x${first_since_master}" = "x" ] ; then
        first_since_master="master"
    fi
    local last_master=${first_since_master}
    local last_master_for_file=${last_master}
    for commit in $(git log --pretty=format:'%h' ${last_master}) ; do
        if git diff --name-only ${commit} ${last_master} | \
                grep -q "${grep_pattern}" ; then
            break
        else
            last_master_for_file=${commit}
        fi
    done
    if [ "x${last_master}" = "x" ] ; then
        last_master_for_file="master"
    fi
    echo ${last_master_for_file}
}


get_branch_name() {
    local branch=$(git rev-parse --abbrev-ref HEAD)
    echo ${branch}
}


write_version_string() {
    local macro=${1}
    local file=${2}
    local version_str=${3}
    echo -e "#ifndef ${macro}\n\
#define ${macro} \"${version_str}\"\n\
#endif" > ${file}
    if [ "${file}" = "src/colvars_version.h" ] ; then
        local version_str_tex=$(version_str_for_tex ${version_str})
        echo -E "\newcommand{\cvversion}{${version_str_tex}}" \
             > doc/cv_version.tex
    fi
}


version_str_for_tex() {
    echo "${1}" | sed 's/_/\\_/'
}


write_version_branch() {
    local function_name=${1}
    local macro=${2}
    local file=${3}
    local branch=${4}
    local version_str=$(${function_name})
    version_str=${version_str%%_*}
    write_version_string ${macro} ${file} "${version_str}_${branch}" 
}


update_version_string() {
    local name=${1}
    local grep_pattern=${2}
    local function_name=${3}
    local macro=${4}
    local file=${5}
    local last_commit=${6}
    
    local branch=$(get_branch_name)
    local version_str=$(date +'%Y-%m-%d')

    if [ "x${last_commit}" = "x" ] ; then
        if [ "${branch}" = "master" ] ; then
            # Just use previous commit
            last_commit=$(git rev-list HEAD^ |head -1)
        else
            # Get the version string from the last commit on master
            last_commit=$(get_last_master_commit "${grep_pattern}")
        fi
    fi
    
    version_str=$(${function_name} ${last_commit})
    # Bump up version when files were modified
    if git diff --name-only ${last_commit} . | \
        grep -q "${grep_pattern}" ; then
        version_str=$(date +'%Y-%m-%d')
    fi
    echo "Setting ${name} version string to ${version_str}"
    write_version_string ${macro} ${file} ${version_str}
    git add ${file}
    git add doc/cv_version.tex
}

get_all_versions() {
    local commit="${1}"
    get_colvarmodule_version "${commit}"
    get_colvarproxy_lammps_version "${commit}"
    get_colvarproxy_namd_version "${commit}"
    get_colvarproxy_vmd_version "${commit}"
}


update_all_versions() {
    local branch="${1}"
    update_version_string "Colvars" \
                          '^src/colvar' \
                          get_colvarmodule_version \
                          COLVARS_VERSION \
                          src/colvars_version.h \
                          ${branch}
    update_version_string "LAMMPS interface" \
                          '^lammps/src/USER-COLVARS/colvarproxy\|^lammps/src/USER-COLVARS/fix_colvars\|^src/colvarproxy' \
                          get_colvarproxy_lammps_version \
                          COLVARPROXY_VERSION \
                          lammps/src/USER-COLVARS/colvarproxy_lammps_version.h \
                          ${branch}
    update_version_string "NAMD interface" \
                          '^namd/src/colvarproxy\|^src/colvarproxy' \
                          get_colvarproxy_namd_version \
                          COLVARPROXY_VERSION \
                          namd/src/colvarproxy_namd_version.h \
                          ${branch}
    update_version_string "VMD interface" \
                          '^vmd/src/colvarproxy\|^src/colvarproxy' \
                          get_colvarproxy_vmd_version \
                          COLVARPROXY_VERSION \
                          vmd/src/colvarproxy_vmd_version.h \
                          ${branch}
}


write_all_versions_branch() {
    local branch="${1}"
    if [ "x${1}" = "x" ] ; then
        echo "Error: write_all_versions_branch requires argument."
        return 1
    fi
    write_version_branch get_colvarmodule_version \
                         COLVARS_VERSION \
                         src/colvars_version.h \
                         ${branch}
    write_version_branch get_colvarproxy_lammps_version \
                         COLVARPROXY_VERSION \
                         lammps/src/USER-COLVARS/colvarproxy_lammps_version.h \
                         ${branch}
    write_version_branch get_colvarproxy_namd_version \
                         COLVARPROXY_VERSION \
                         namd/src/colvarproxy_namd_version.h \
                         ${branch}
    write_version_branch get_colvarproxy_vmd_version \
                         COLVARPROXY_VERSION \
                         vmd/src/colvarproxy_vmd_version.h \
                         ${branch}
}
