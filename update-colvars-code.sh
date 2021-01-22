#!/bin/bash
# -*- sh-basic-offset: 2; sh-indentation: 2; -*-

# Script to update a NAMD, VMD, LAMMPS or GROMACS source tree with the latest Colvars
# version.

# Enforce using portable C locale
LC_ALL=C
export LC_ALL

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

# Was the target Makefile changed?
updated_makefile=0

# Was the last file updated?
updated_file=0

force_update=0
if [ $1 = "-f" ]
then
  echo "Forcing update of all files"
  force_update=1
  shift
fi

# Undocumented flag
reverse=0
if [ $1 = "-R" ]
then
  echo "Reverse: updating git tree from downstream tree"
  reverse=1
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
elif [ -f "${target}/src/gromacs/commandline.h" ]
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
  elif [ -f "${target}/src/gromacs/commandline.h" ]
  then
    code="GROMACS"
  else
    echo "ERROR: Cannot detect a supported code in the target directory."
    exit 3
  fi
fi

get_gromacs_major_version_cmake() {
  cat $1 | grep 'set(GMX_VERSION_MAJOR' | \
    sed -e 's/set(GMX_VERSION_MAJOR //' -e 's/)//'
}

get_gromacs_minor_version_cmake() {
  cat $1 | grep 'set(GMX_VERSION_PATCH' | \
    sed -e 's/set(GMX_VERSION_PATCH //' -e 's/)//'
}

echo "Detected ${code} source tree in ${target}"
if [ ${code} = "GROMACS" ]
then
  GMX_VERSION_INFO=${target}/cmake/gmxVersionInfo.cmake
  if [ ! -f ${GMX_VERSION_INFO} ] ; then
    echo "ERROR: Cannot find file ${GMX_VERSION_INFO}."
    exit 3
  fi

  GMX_MAJOR_VERSION=`get_gromacs_major_version_cmake ${GMX_VERSION_INFO}`
  GMX_MINOR_VERSION=`get_gromacs_minor_version_cmake ${GMX_VERSION_INFO}`

  GMX_VERSION=${GMX_MAJOR_VERSION}.${GMX_MINOR_VERSION}
  echo "Detected GROMACS version ${GMX_VERSION}."

  case ${GMX_VERSION} in
    2018*)
      GMX_VERSION='2018.x'
      ;;
    2020*)
      GMX_VERSION='2020.x'
      ;;
    *)
    GMX_VERSION="master"
    if [ $force_update = 0 ] ; then
      echo " ******************************************************************************"
      echo "  ERROR: Support for GROMACS version ${GMX_VERSION} has not been tested."
      echo "  You may override with -f, but be mindful of compilation or runtime problems."
      echo " ******************************************************************************"
      exit 3
    fi
    ;;
  esac

fi
echo -n "Updating ..."


# Conditional file copy
condcopy() {
  if [ $reverse -eq 1 ]
  then
    a=$2
    b=$1
    PATCH_OPT="-R"
  else
    a=$1
    b=$2
    PATCH_OPT=""
  fi

  updated_file=0

  if [ -d $(dirname "$b") ]
  then
    if [ $checkonly -eq 1 ]
    then
      cmp -s "$a" "$b" || diff -uNw "$b" "$a"
    else
      if ! cmp -s "$a" "$b" ; then
        cp "$a" "$b"
        updated_file=1
      fi
      echo -n '.'
    fi
  fi
}


# Check files related to, but not part of the Colvars module
checkfile() {
  if [ $reverse -eq 1 ]
  then
    a=$2
    b=$1
  else
    a=$1
    b=$2
  fi
  diff -uNw "${a}" "${b}" > $(basename ${a}).diff
  if [ -s $(basename ${a}).diff ]
  then
    echo "Differences found between ${a} and ${b} -- Check $(basename ${a}).diff and merge changes as needed, or use the -f flag."
    if [ $force_update = 1 ]
    then
      echo "Overwriting ${b}, as requested by the -f flag."
      cp "$a" "$b"
    fi
  else
    rm -f $(basename ${a}).diff
  fi
}


# Update LAMMPS tree
if [ ${code} = "LAMMPS" ]
then

  # Update code-independent headers and sources
  for src in ${source}/src/colvar*.h ${source}/src/colvar*.cpp
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/lib/colvars/${tgt}"
  done

  # Update makefiles for library
  for src in ${source}/lammps/lib/colvars/Makefile.{common,deps}
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/lib/colvars/${tgt}"
  done

  # Update LAMMPS interface files (package part)
  if [ -f ${target}/src/random_park.h ]
  then
    # GitHub version, using old pseudo random number generators
    for src in ${source}/lammps/src/USER-COLVARS/colvarproxy_lammps.cpp \
               ${source}/lammps/src/USER-COLVARS/colvarproxy_lammps.h \
               ${source}/lammps/src/USER-COLVARS/colvarproxy_lammps_version.h \
               ${source}/lammps/src/USER-COLVARS/fix_colvars.cpp \
               ${source}/lammps/src/USER-COLVARS/fix_colvars.h
    do \
      tgt=$(basename ${src})
      condcopy "${src}" "${target}/src/USER-COLVARS/${tgt}"
    done
  else
    echo "ERROR: Support for the new pRNG (old LAMMPS-ICMS branch) is currently disabled."
    exit 2
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
    condcopy "${src}" "${docdir}/src/PDF/${tgt}"
  done

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

  # New layout: copy library files to the "colvars" folder
  for src in ${source}/src/*.h ${source}/src/*.cpp
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/colvars/src/${tgt}"
  done
  condcopy "${source}/namd/colvars/src/Makefile.namd" \
           "${target}/colvars/src/Makefile.namd"
  if [ $updated_file = 1 ] ; then
    updated_makefile=1
  fi
  condcopy "${source}/namd/colvars/Make.depends" \
           "${target}/colvars/Make.depends"

  # Update abf_integrate
  for src in ${source}/colvartools/*h ${source}/colvartools/*cpp
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/lib/abf_integrate/${tgt}"
  done
  condcopy "${source}/colvartools/Makefile" \
           "${target}/lib/abf_integrate/Makefile"

  # Update NAMD interface files
  for src in \
      ${source}/namd/src/colvarproxy_namd.h \
      ${source}/namd/src/colvarproxy_namd_version.h \
      ${source}/namd/src/colvarproxy_namd.C
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/src/${tgt}"
  done

  # Is this a devel branch of NAMD 3?
  if echo $NAMD_VERSION | grep -q '3.0a'
  then
    echo "Detected a devel version of NAMD 3:"
    echo "Assuming version number below 2.14b1 to disable Volmaps"
    sed -i 's/\#define\ NAMD_VERSION_NUMBER\ 34471681/\#define\ NAMD_VERSION_NUMBER\ 34471680/' ${target}/src/colvarproxy_namd.h
  fi

  # Update replacement text for the Colvars manual
  condcopy "${source}/namd/ug/ug_colvars.tex" \
           "${target}/ug/ug_colvars.tex"

  echo ' done.'

  # Check for changes in related NAMD files
  for src in \
      ${source}/namd/src/GlobalMasterColvars.h \
      ${source}/namd/src/ScriptTcl.h \
      ${source}/namd/src/ScriptTcl.C \
      ${source}/namd/src/SimParameters.h \
      ${source}/namd/src/SimParameters.C \
      ;
  do \
    tgt=$(basename ${src})
    checkfile "${src}" "${target}/src/${tgt}"
  done
  for src in ${source}/namd/Make* ${source}/namd/config
  do
    tgt=$(basename ${src})
    checkfile "${src}" "${target}/${tgt}"
    if [ $updated_file = 1 ] ; then
      updated_makefile=1
    fi
  done

  # One last check that each file is correctly included in the dependencies
  for file in ${target}/colvars/src/*.{cpp,h} ; do
    if [ ! -f ${target}/colvars/Make.depends ] || \
       [ ! -f ${target}/lepton/Make.depends ] ; then
      updated_makefile=1
      break
    fi
    if ! grep -q ${file} ${target}/colvars/Make.depends ; then
      updated_makefile=1
    fi
  done

  if [ $updated_makefile = 1 ] ; then
    echo ""
    echo "  *************************************************"
    echo "    Please run \"make depends\" in the NAMD tree."
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

  condcopy "${source}/vmd/src/colvars_files.pl" "${target}/src/colvars_files.pl"

  echo ' done.'

  exit 0
fi


# Update VMD plugins tree
if [ ${code} = "VMD-PLUGINS" ]
then

  # Use the Dashboard's Makefile to patch the plugin tree
  if pushd ${source}/vmd/cv_dashboard > /dev/null ; then
    DASHBOARD_VERSION=$(grep ^VERSION Makefile.local | cut -d' ' -f 3)
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

# Update GROMACS tree
if [ ${code} = "GROMACS" ]
then
  echo ""
  if [ -d ${target}/src/gromacs/colvars ]
  then
    echo "Your ${target} source tree seems to have already been patched."
    echo "Update with the last Colvars source."
  else
    mkdir ${target}/src/gromacs/colvars
  fi

  # copy library files and proxy files to the "colvars" folder
  for src in ${source}/src/*.h ${source}/src/*.cpp ${source}/gromacs/src/*.h ${source}/gromacs/gromacs-${GMX_VERSION}/*
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/src/gromacs/colvars/${tgt}" "${cpp_patch}"
  done
  echo ""

  # Apply patch for Gromacs files
  patch --forward -s -p0 -d ${target} < ${source}/gromacs/gromacs-${GMX_VERSION}.patch
  ret_val=$?
  if [ $ret_val -ne 0 ]
  then
    echo " ************************************************************************* "
    echo " Patch fails. It seems the Gromacs source files have been already patched. "
    echo " ************************************************************************* "
  else
    echo ' done.'
    echo ""
    echo "  *******************************************"
    echo "    Please create your build with cmake now."
    echo "  *******************************************"
  fi

  # Update the proxy version if needed
  shared_gmx_proxy_version=$(grep '^#define' "${source}/gromacs/src/colvarproxy_gromacs_version.h" | cut -d' ' -f 3)

  patch_gmx_proxy_version=$(grep '^#define' "${target}/src/gromacs/colvars/colvarproxy_gromacs_version.h" | cut -d' ' -f 3)

  if [ ${shared_gmx_proxy_version} \> ${patch_gmx_proxy_version} ] ; then
    condcopy ${source}/gromacs/src/colvarproxy_gromacs_version.h \
      "${target}/src/gromacs/colvars/colvarproxy_gromacs_version.h"
  fi

  exit 0
fi
