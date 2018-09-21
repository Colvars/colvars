#!/bin/sh
# -*- sh-basic-offset: 2; sh-indentation: 2; -*-

# Script to update a NAMD, VMD, or LAMMPS source tree with the latest Colvars
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
   supported MD codes: NAMD, VMD, LAMMPS

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

cpp_patch=${source}/devel-tools/update-header-cpp.patch
tex_patch=${source}/devel-tools/update-header-tex.patch

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
    echo ERROR: Cannot detect a supported code in the target directory
    exit 3
  fi
fi

echo "Detected ${code} source tree in ${target}"
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

  TMPFILE=`mktemp`

  # If a patch file is available, apply it to the source file
  # (reversed if necessary)
  if [ "x$3" != "x" ] ; then
    if [ -f "$3" ] ; then
      patch $PATCH_OPT < $3 $a -o $TMPFILE > /dev/null
      # Patched file is new source
      a=$TMPFILE
    fi
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

  rm -f $TMPFILE
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
  for src in ${source}/src/*.h ${source}/src/*.cpp
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/lib/colvars/${tgt}" "${cpp_patch}"
  done

  # Update makefiles for library
  for src in \
    ${source}/lammps/lib/colvars/Makefile.* \
    ${source}/lammps/lib/colvars/Install.py \
    ${source}/lammps/lib/colvars/README
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
      condcopy "${src}" "${target}/src/USER-COLVARS/${tgt}" "${cpp_patch}"
    done
    for src in ${source}/lammps/src/USER-COLVARS/Install.sh \
               ${source}/lammps/src/USER-COLVARS/group_ndx.cpp \
               ${source}/lammps/src/USER-COLVARS/group_ndx.h \
               ${source}/lammps/src/USER-COLVARS/ndx_group.cpp \
               ${source}/lammps/src/USER-COLVARS/ndx_group.h \
               ${source}/lammps/src/USER-COLVARS/README
    do \
      tgt=$(basename ${src})
      condcopy "${src}" "${target}/src/USER-COLVARS/${tgt}"
    done
  else
    echo "ERROR: Support for the new pRNG (old LAMMPS-ICMS branch) is currently disabled."
    exit 2
  fi

  # Update LAMMPS documentation
  # Location of documentation has changed with version 10 May 2016
  test -d "${target}/doc/src/PDF" && \
    docdir="${target}/doc/src" || docdir="${target}/doc"
  for src in ${source}/lammps/doc/*.txt
    do \
      tgt=$(basename ${src})
    condcopy "${src}" "${docdir}/${tgt}"
  done

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
    condcopy "${src}" "${docdir}/PDF/${tgt}"
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

  if [ ! -d ${target}/colvars ] ; then
    # Old layout: copy library files to the "src" folder
    echo ""
    if [ $force_update = 1 ] ; then
      echo "Creating folder containing Colvars library files"
      mkdir ${target}/colvars
      mkdir ${target}/colvars/src
      rm -f ${target}/src/colvar*
    else
      echo "Error: Colvars files are now stored in a separate folder."
      echo "Please upgrade to the latest NAMD version (supported), "
      echo "or use the -f flag (unsupported)."
      exit 1
    fi
  fi
  
  # New layout: copy library files to the "colvars" folder
  for src in ${source}/src/*.h ${source}/src/*.cpp
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/colvars/src/${tgt}" "${cpp_patch}"
  done
  condcopy "${source}/namd/colvars/src/Makefile.namd" \
           "${target}/colvars/src/Makefile.namd"
  if [ $updated_file = 1 ] ; then
    updated_makefile=1
  fi

  # Update abf_integrate
  for src in ${source}/colvartools/*h ${source}/colvartools/*cpp
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/lib/abf_integrate/${tgt}" "${cpp_patch}"
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
    condcopy "${src}" "${target}/src/${tgt}" "${cpp_patch}"
  done

  # Copy doc files
  condcopy "${source}/doc/colvars-refman.bib" \
           "${target}/ug/ug_colvars.bib"
  condcopy "${source}/doc/colvars-refman-main.tex" \
           "${target}/ug/ug_colvars.tex" "${tex_patch}"
  condcopy "${source}/doc/colvars-cv.tex" \
           "${target}/ug/ug_colvars-cv.tex" "${tex_patch}"
  condcopy "${source}/namd/ug/ug_colvars_macros.tex" \
           "${target}/ug/ug_colvars_macros.tex" "${tex_patch}"
  condcopy "${source}/doc/colvars_diagram.pdf" \
           "${target}/ug/figures/colvars_diagram.pdf"
  condcopy "${source}/doc/colvars_diagram.eps" \
           "${target}/ug/figures/colvars_diagram.eps"

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
    condcopy "${src}" "${target}/src/${tgt}" "${cpp_patch}"
  done
  # Update code-independent sources
  for src in ${source}/src/*.cpp
  do \
    tgt=$(basename ${src%.cpp})
    condcopy "${src}" "${target}/src/${tgt}.C" "${cpp_patch}"
  done

  condcopy "${source}/doc/colvars-refman.bib" \
           "${target}/doc/ug_colvars.bib"
  condcopy "${source}/doc/colvars-refman-main.tex" \
           "${target}/doc/ug_colvars.tex" "${tex_patch}"
  condcopy "${source}/doc/colvars-cv.tex" \
           "${target}/doc/ug_colvars-cv.tex" "${tex_patch}"
  condcopy "${source}/vmd/doc/ug_colvars_macros.tex" \
           "${target}/doc/ug_colvars_macros.tex" "${tex_patch}"
  condcopy "${source}/doc/colvars_diagram.pdf" \
           "${target}/doc/pictures/colvars_diagram.pdf"
  condcopy "${source}/doc/colvars_diagram.eps" \
           "${target}/doc/pictures/colvars_diagram.eps"

  # Update VMD interface files
  for src in \
      ${source}/vmd/src/colvarproxy_vmd.h \
      ${source}/vmd/src/colvarproxy_vmd_version.h \
      ${source}/vmd/src/colvarproxy_vmd.C  
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/src/${tgt}" "${cpp_patch}"
  done

  condcopy "${source}/vmd/src/colvars_files.pl" "${target}/src/colvars_files.pl"

  echo ' done.'

  # Check for changes in related VMD files
  for src in ${source}/vmd/src/tcl_commands.C
  do \
    tgt=$(basename ${src})
    checkfile "${src}" "${target}/src/${tgt}"
  done
  for src in ${source}/vmd/configure
  do 
    tgt=$(basename ${src})
    checkfile "${src}" "${target}/${tgt}"
  done

  exit 0
fi


# Update GROMACS tree
if [ ${code} = "GROMACS" ]
then
  echo "Error: the GROMACS implementation of Colvars is not under active development."
  exit 1
fi
