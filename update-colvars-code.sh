#!/bin/bash
# -*- sh-basic-offset: 2; sh-indentation: 2; -*-

set -e

# Script to update a NAMD, VMD, LAMMPS or GROMACS source tree with the latest Colvars
# version.

# Enforce using portable C locale
LC_ALL=C
export LC_ALL

# Detect Git
if [ -z "${GIT}" ] && hash git >& /dev/null ; then
  GIT=$(hash -t git)
fi
if [ -n "${GIT}" ] ; then
  if [ ! -x "${GIT}" ] ; then
    echo "Error: \$GIT variable defined, but its value is not an executable: \"$GIT\"" >&2
    exit 1
  fi
fi

if [ $# -lt 1 ]
then
    cat <<EOF

 usage: sh $0 [-f] <target source tree>

   -f  "force-update": overwrite conflicting files such as Makefile
        (default: create diff files for inspection --- MD code may be different)

   <target source tree> = root directory of the MD code sources
   supported codes: NAMD, VMD, VMD PLUGINS, LAMMPS, GROMACS

EOF
   exit 1
fi


force_update=0
if [ $1 = "-f" ]
then
  echo "Forcing update of all files"
  force_update=1
  shift
fi

# Assume forward patching
PATCH_OPTS="-p1 --forward -s -r -"

# Undocumented flag
reverse=0
if [ $1 = "-R" ]
then
  echo "Reverse: updating git tree from downstream tree"
  reverse=1
  PATCH_OPTS="-p1 --reverse -s"
  shift
fi

# Infer source path from name of script
source=$(dirname "$0")

# Check general validity of target path
target="$1"
if [ ! -d "${target}" ]
then
    echo "ERROR: Target directory ${target} does not exist"
    exit 2
fi

# Undocumented option to only compare trees
checkonly=0
[ "$2" = "--diff" ] && checkonly=1
[ $force_update = 1 ] && checkonly=0

# Try to determine what code resides inside the target dir
code=unknown
if [ -f "${target}/src/lammps.h" ]
then
  code="LAMMPS"
elif [ -f "${target}/src/NamdTypes.h" ]
then
  code="NAMD"
elif [ -f "${target}/src/VMDApp.h" ]
then
  code="VMD"
elif [ -f "${target}/include/molfile_plugin.h" ]
then
  code="VMD-PLUGINS"
elif [ -f "${target}/src/gromacs/version.h.cmakein" ]
then
  code="GROMACS"
else
  # Handle the case if the user points to ${target}/src
  target=$(dirname "${target}")
  if [ -f "${target}/src/lammps.h" ]
  then
    code="LAMMPS"
  elif [ -f "${target}/src/NamdTypes.h" ]
  then
    code="NAMD"
  elif [ -f "${target}/src/VMDApp.h" ]
  then
    code="VMD"
  elif [ -f "${target}/src/gromacs/version.h.cmakein" ]
  then
    code="GROMACS"
  else
    echo "ERROR: Cannot detect a supported code in the target directory."
    exit 3
  fi
fi

if [ ${code} == "GROMACS" ] ; then
  if [ ! -f ${target}/src/gromacs/applied_forces/colvars/colvarsMDModule.cpp ] ; then
    echo "ERROR: Versions of GROMACS prior to 2024 are no longer supported."
    exit 1
  fi
fi

COLVARS_VERSION=$(grep '^#define COLVARS_VERSION' $(dirname $0)/src/colvars_version.h | cut -d' ' -f 3 | tr -d '"')
if [ -z "${COLVARS_VERSION}" ] ; then
  echo "Error reading Colvars version." >&2
  exit 1
fi


get_gromacs_major_version_cmake() {
  cat $1 | grep 'set(GMX_VERSION_MAJOR' | \
    sed -e 's/set(GMX_VERSION_MAJOR //' -e 's/)//'
}

get_gromacs_minor_version_cmake() {
  cat $1 | grep 'set(GMX_VERSION_PATCH' | \
    sed -e 's/set(GMX_VERSION_PATCH //' -e 's/)//'
}


UPDATE_LEPTON=${UPDATE_LEPTON:-yes}

copy_lepton() {

  local target_path=${1}

  if [ -d ${source}/openmm-source ] ; then
    OPENMM_SOURCE=${source}/openmm-source
  elif [ -d ${source}/../openmm-source ] ; then
    OPENMM_SOURCE=${source}/../openmm-source
  fi

  if [ -z "${OPENMM_SOURCE}" ] ; then
    # Create a temp folder and download OpenMM into it
    OPENMM_SOURCE=$(mktemp -d /tmp/openmm-source-XXXXXX)
  fi

  # Download Lepton if needed
  if [ ! -d ${OPENMM_SOURCE}/libraries/lepton ] ; then
    if [ -n "${GIT}" ] ; then
      echo "Downloading Lepton library (used in Colvars) via the OpenMM repository"
      ${GIT} clone --depth=1 https://github.com/openmm/openmm.git ${OPENMM_SOURCE}
    fi
  fi

  # Copy Lepton into target source tree
  if [ -d ${OPENMM_SOURCE}/libraries/lepton ] ; then
    cp -f -p -R ${OPENMM_SOURCE}/libraries/lepton ${target_path}
  else
    echo "ERROR: could not download the Lepton library automatically." >&2
    echo "       Please clone the OpenMM repository (https://github.com/openmm/openmm) " >&2
    echo "       in a directory of your choice, and set the environment variable OPENMM_SOURCE " >&2
    echo "       to the absolute path of that directory." >&2
    return 1
  fi
}


echo "Detected ${code} source tree in ${target}"
if [ ${code} = "GROMACS" ]
then
  GMX_VERSION_INFO=${target}/cmake/gmxVersionInfo.cmake
  if [ ! -f ${GMX_VERSION_INFO} ] ; then
    echo "ERROR: Cannot find file ${GMX_VERSION_INFO}."
    exit 3
  fi

  CMAKE_LISTS=${target}/CMakeLists.txt
  if [ ! -f ${CMAKE_LISTS} ] ; then
    echo "ERROR: Cannot find file ${CMAKE_LISTS}."
    exit 3
  fi
  GMX_VERSION=$(cat ${CMAKE_LISTS} | grep '^project(' | sed -e 's/project(Gromacs VERSION //' -e 's/)//')

  echo "Detected GROMACS version ${GMX_VERSION}."

  if [ -z "${GITHUB_ACTION}" ] ; then
    # Only set version outside CI to avoid invalidating the compiler cache
    if grep -q 'set(GMX_VERSION_STRING_OF_FORK ""' ${GMX_VERSION_INFO} ; then
      sed -i "s/set(GMX_VERSION_STRING_OF_FORK \"\"/set(GMX_VERSION_STRING_OF_FORK \"Colvars-${COLVARS_VERSION}\"/" ${GMX_VERSION_INFO}
    fi
  fi

fi

echo -n "Updating ..."


# Conditional file copy
condcopy() {
  if [ $reverse -eq 1 ]
  then
    a=$2
    b=$1
  else
    a=$1
    b=$2
  fi

  if [ -d $(dirname "$b") ]
  then
    if [ $checkonly -eq 1 ]
    then
      cmp -s "$a" "$b" || diff -uNw "$b" "$a"
    else
      if ! cmp -s "$a" "$b" ; then
        cp "$a" "$b"
      fi
      echo -n '.'
    fi
  fi
}


# Set updated_makefile to 1 if file is changed
check_target_makefile() {
  local target=$1 file=$2
  if [ -z "${updated_makefile}" ] ; then
    updated_makefile=0
  fi
  if [ -n "${GIT}" ] ; then
    if ! git -C ${target} diff --quiet ${file} ; then
      echo -n ".(note: file ${file} was modified)"
      updated_makefile=$((${updated_makefile} || 1))
    fi
  else
    updated_makefile=1
  fi
}


# Update LAMMPS tree
if [ ${code} = "LAMMPS" ]
then

  # # LAMMPS contains an updated but patched Lepton, uncomment here if needed
  # copy_lepton ${target}/lib/ || exit 1

  # Update code-independent headers and sources
  for src in ${source}/src/colvar*.h ${source}/src/colvar*.cpp
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/lib/colvars/${tgt}"
  done
  # Update the Colvars CUDA source files
  mkdir -p "${target}/lib/colvars/cuda"
  for src in ${source}/src/cuda/*.h ${source}/src/cuda/*.cu
  do
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/lib/colvars/cuda/${tgt}"
  done

  for src in \
    ${source}/lammps/src/COLVARS/colvarproxy_lammps{.cpp,.h,_version.h} \
    ${source}/lammps/src/COLVARS/fix_colvars.{cpp,h} \
    ${source}/lammps/src/COLVARS/inthash.{cpp,h}
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/src/COLVARS/${tgt}"
  done

  if [ -f ${source}/lammps/COLVARS.cmake.patch ] ; then
    # Do not exit if the patch fails - already applied in development branch
    patch ${PATCH_OPTS} -d ${target} < ${source}/lammps/COLVARS.cmake.patch || true
  fi

  downloaded_pdf=0
  # Copy PDF of the user manual
  if [ ! -f ${source}/doc/colvars-refman-lammps.pdf ] ; then
    if curl -L -o ${source}/doc/colvars-refman-lammps.pdf \
            https://colvars.github.io/pdf/colvars-refman-lammps.pdf \
        1> /dev/null 2> /dev/null || \
        wget -O ${source}/doc/colvars-refman-lammps.pdf \
              https://colvars.github.io/pdf/colvars-refman-lammps.pdf \
        1> /dev/null 2> /dev/null \
       ; then
      downloaded_pdf=1
      echo -n '.'
    else
      echo ""
      echo "Error: could not download the PDF manual automatically."
      echo "Please download it manually from:"
      echo "  https://colvars.github.io/pdf/colvars-refman-lammps.pdf"
      echo "and copy it into ${source}/doc,"
      echo "or re-generate it using:"
      echo "  cd ${source}/doc ; make colvars-refman-lammps.pdf; cd -"
      exit 1
    fi
  fi
  for src in ${source}/doc/colvars-refman-lammps.pdf
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/doc/src/PDF/${tgt}"
  done
  condcopy "${source}/lammps/doc/src/fix_colvars.rst" \
           "${target}/doc/src/fix_colvars.rst"
  rm -f "${target}/doc/src/fix_colvars.txt"

  echo ' done.'
  if [ ${downloaded_pdf} = 1 ] ; then
    echo "Note: the PDF manual for the latest Colvars version was downloaded.  "
    echo "If you are using an older version, you can generate the corresponding PDF with:"
    echo "  cd ${source}/doc ; make colvars-refman-lammps.pdf; cd -"
    echo "and run this script a second time."
  fi
  exit 0
fi


# Update NAMD tree
if [ ${code} = "NAMD" ]
then
  NAMD_VERSION=$(grep ^NAMD_VERSION ${target}/Makefile | cut -d' ' -f3)

  updated_makefile=0

  if [ "x${UPDATE_LEPTON}" == "xyes" ] ; then
    echo -n "(note: updating Lepton)"

    if ! copy_lepton ${target}/ ; then
      echo "Error: in updating Lepton" >&2
      exit 1
    fi

    condcopy "${source}/namd/lepton/Make.depends" "${target}/lepton/Make.depends"
    check_target_makefile "${target}" "lepton/Make.depends"

    condcopy "${source}/namd/lepton/Makefile.namd" "${target}/lepton/Makefile.namd"
    check_target_makefile "${target}" "lepton/Makefile.namd"
  fi

  # Copy library files to the "colvars" folder
  for src in ${source}/src/*.h ${source}/src/*.cpp
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/colvars/src/${tgt}"
  done

  # Generate Makefile for library files
  ${source}/devel-tools/generate-namd-makefile.sh \
           COLVARSSRCS '$(COLVARSSRCDIR)' \
           COLVARSLIB '$(DSTDIR)' \
           ${source}/src/*.cpp \
           > "${target}/colvars/src/Makefile.namd"
  check_target_makefile "${target}" "colvars/src/Makefile.namd"

  # Optionally copy over the Colvars library's Make.depends
  if [ -s "${source}/namd/colvars/Make.depends" ] ; then
    condcopy "${source}/namd/colvars/Make.depends" "${target}/colvars/Make.depends"
    check_target_makefile "${target}" "colvars/Make.depends"
  fi

  # Update NAMD interface files
  for src in \
      ${source}/namd/src/*.h \
      ${source}/namd/src/*.C
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/src/${tgt}"
  done
  for patch in ${source}/namd/*.patch ${source}/namd/src/*.patch ; do
    if [ -s ${patch} ] ; then
      patch ${PATCH_OPTS} -d "${target}" < "${patch}" || true
    fi
  done

  # Update the Colvars CUDA source files
  # echo "Copy Colvars CUDA source files"
  mkdir -p "${target}/colvars/src/cuda"
  for src in ${source}/src/cuda/*.h ${source}/src/cuda/*.cu
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/colvars/src/cuda/${tgt}"
  done

  # Update the NAMD CudaGlobalMaster interface files
  # echo "Copy Colvars-NAMD CudaGlobalMaster interface files"
  mkdir -p "${target}/colvars/namd/cudaglobalmaster"
  condcopy "${source}/namd/cudaglobalmaster/Make.depends" \
           "${target}/colvars/namd/cudaglobalmaster/Make.depends"
  condcopy "${source}/namd/cudaglobalmaster/Makefile_common.namd" \
           "${target}/colvars/namd/cudaglobalmaster/Makefile_common.namd"
  condcopy "${source}/namd/cudaglobalmaster/Makefile_cuda.namd" \
           "${target}/colvars/namd/cudaglobalmaster/Makefile_cuda.namd"
  for src in \
      ${source}/namd/cudaglobalmaster/*.h \
      ${source}/namd/cudaglobalmaster/*.C \
      ${source}/namd/cudaglobalmaster/*.cu \
      ${source}/namd/cudaglobalmaster/CMakeLists.txt
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/colvars/namd/cudaglobalmaster/${tgt}"
  done

  # Update abf_integrate
  for src in ${source}/colvartools/*h ${source}/colvartools/*cpp
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/lib/abf_integrate/${tgt}"
  done
  condcopy "${source}/colvartools/Makefile" \
           "${target}/lib/abf_integrate/Makefile"

  # Update replacement text for the Colvars manual
  condcopy "${source}/namd/ug/ug_colvars.tex" \
           "${target}/ug/ug_colvars.tex"

  echo ' done.'

  # One last check that each file is correctly included in the dependencies
  for file in ${target}/colvars/src/*.{cpp,h} ; do
    if ! grep -q ${file#${target}/} ${target}/colvars/Make.depends ; then
      echo "File ${file#${target}/} is missing from the Colvars library's Make.depends"
      updated_makefile=1
    fi
  done

  if [ $updated_makefile != 0 ] ; then
    echo ""
    echo "  *************************************************"
    if [ -n "${GIT}" ] ; then
      echo "    NAMD Makefiles were modified from their Git version:"
    else
      echo "    NAMD Makefiles were probably modified:"
    fi
    echo "    please run \"make depends\" in the NAMD tree."
    echo "  *************************************************"
  fi

  exit 0
fi


# Update VMD tree
if [ ${code} = "VMD" ]
then

  # Update code-independent headers
  for src in ${source}/src/*.h
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/src/${tgt}"
  done
  # Update code-independent sources
  for src in ${source}/src/*.cpp
  do \
    tgt=$(basename ${src%.cpp})
    condcopy "${src}" "${target}/src/${tgt}.C"
  done
  # Update the Colvars CUDA code
  mkdir -p "${target}/src/cuda"
  for src in ${source}/src/cuda/*.h ${source}/src/cuda/*.cu
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/src/cuda/${tgt}"
  done

  if [ -s ${source}/vmd/src/tcl_commands.C.patch ] ; then
    # Do not exit if the patch fails - already applied in development branch
    patch ${PATCH_OPTS} -d ${target} < ${source}/vmd/src/tcl_commands.C.patch || true
  fi

  # Update replacement text for the Colvars manual
  condcopy "${source}/vmd/doc/ug_colvars.tex" \
           "${target}/doc/ug_colvars.tex"

  # Update VMD interface files
  for src in \
      ${source}/vmd/src/colvarproxy_vmd.h \
      ${source}/vmd/src/colvarproxy_vmd_version.h \
      ${source}/vmd/src/colvarproxy_vmd.C
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/src/${tgt}"
  done

  bash "${source}/devel-tools/generate-vmd-makefile.sh" \
       "${source}/src/"*.cpp \
       "${source}/src/"*.h \
       "${source}/vmd/src/"*.C \
       "${source}/vmd/src/"*.h \
       > "${target}/src/colvars_files.pl"

  echo ' done.'

  exit 0
fi


# Update VMD plugins tree
if [ ${code} = "VMD-PLUGINS" ]
then

  # Use the Dashboard's Makefile to patch the plugin tree
  if pushd ${source}/vmd/cv_dashboard > /dev/null ; then
    DASHBOARD_VERSION=$(cat VERSION)
    if [ -d ${target}/noarch ] ; then
      # This is an already-installed plugin tree
      DASHBOARD_DESTINATION=${target}/noarch/tcl/cv_dashboard${DASHBOARD_VERSION}
    else
      # This is the source tree
      DASHBOARD_DESTINATION=${target}/cv_dashboard
    fi
    DESTINATION=${DASHBOARD_DESTINATION} \
      make --quiet -f Makefile.local > /dev/null
    echo -n '......'
    popd > /dev/null
  fi
  echo ' done.'
fi


# Update GROMACS tree (MDModules interface-ready version)
if [ ${code} = "GROMACS" ]
then

  target_folder=${target}/src/external/colvars

  mkdir -p ${target_folder}

  # Copy library files and proxy files to the "src/external/colvars" folder
  for src in ${source}/src/*.h ${source}/src/*.cpp
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target_folder}/${tgt}"
  done

  # Update the Colvars CUDA source files
  mkdir -p "${target_folder}/cuda"
  for src in ${source}/src/cuda/*.h ${source}/src/cuda/*.cu
  do
    tgt=$(basename ${src})
    condcopy "${src}" "${target_folder}/cuda/${tgt}"
  done
  echo ""

  # Patch CMake build recipes when applicable
  if [ -s ${source}/gromacs/gmxManageColvars.cmake.patch ] ; then
    patch ${PATCH_OPTS} -d ${target} < ${source}/gromacs/gmxManageColvars.cmake.patch || true
  fi
  if [ -s ${source}/gromacs/CMakeLists.txt.patch ] ; then
    patch ${PATCH_OPTS} -d ${target} < ${source}/gromacs/CMakeLists.txt.patch || true
  fi
  if [ -s ${source}/gromacs/CMakeLists.txt.patch ] ; then
    patch ${PATCH_OPTS} -d ${target} < ${source}/gromacs/CMakeLists.txt.patch || true
  fi

  if [ "x${UPDATE_LEPTON}" == "xyes" ] ; then
    echo -n "(note: adding/updating Lepton)"
    copy_lepton ${target}/src/external/ || exit 1
    if [ -s ${source}/gromacs/CMakeLists.txt.Lepton.patch ] ; then
      patch ${PATCH_OPTS} -d ${target} < ${source}/gromacs/CMakeLists.txt.Lepton.patch || true
    fi
    condcopy ${source}/gromacs/cmake/gmxManageLepton.cmake "${target}/cmake/gmxManageLepton.cmake"
  fi

  echo

  # Copy MDModules files to the "src/gromacs/src/applied_forces/colvars" folder
  target_folder=${target}/src/gromacs/applied_forces/colvars

  for patch_file in ${source}/gromacs/src/applied_forces/colvars/*.patch ; do
    if [ -s ${patch_file} ] ; then
      # Do not exit if the patch fails - already applied
      patch ${PATCH_OPTS} -d ${target} < ${patch_file} || true
    fi
  done

  if [ -d ${target_folder} ]
  then
    echo "${target} source tree seems to have already been patched."
    echo "Preserving files from the Gromacs tree."
  else
    mkdir ${target_folder}
    mkdir -p ${target_folder}/tests/refdata
    condcopy gromacs/src/colvarproxy_gromacs_version.h "${target_folder}/colvarproxy_gromacs_version.h"
    for src in ${source}/gromacs/src/applied_forces/colvars/*.* ; do
      tgt=$(basename ${src})
      condcopy "${src}" "${target_folder}/${tgt}"
    done
    for src in ${source}/gromacs/src/applied_forces/colvars/tests/*.* ; do
      tgt=$(basename ${src})
      condcopy "${src}" "${target_folder}/tests/${tgt}"
    done
    for src in ${source}/gromacs/src/applied_forces/colvars/tests/refdata/*.* ; do
      tgt=$(basename ${src})
      condcopy "${src}" "${target_folder}/tests/refdata/${tgt}"
    done
  fi

  echo ""

  # Apply patch for Gromacs files
  if [ -s ${source}/gromacs/gromacs-mdmodules.patch ] ; then
    patch ${PATCH_OPTS} -d ${target} < ${source}/gromacs/gromacs-mdmodules.patch || true
  fi
  ret_val=$?
  if [ $ret_val -ne 0 ]
  then
    echo " ************************************************************************************************ "
    echo " Patch fails. Your GROMACS developement tree may be different of one used for creating the patch. "
    echo " ************************************************************************************************ "
  else
    echo ' done.'
    echo ""
    echo "  *******************************************"
    echo "    Please create your build with cmake now."
    echo "  *******************************************"
  fi

  exit 0
fi
