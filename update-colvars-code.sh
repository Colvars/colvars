#!/bin/sh
# Script to update a NAMD, VMD, LAMMPS, or GROMACS source tree with the latest colvars code.

if [ $# -lt 1 ]
then
    cat <<EOF

 usage: sh $0 [-f] <target source tree>

   -f  "force-update": overwrite conflicting files such as Makefile
        (default: create diff files for inspection --- MD code may be different)

   <target source tree> = root directory of the MD code sources
   supported MD codes: NAMD, VMD, LAMMPS, GROMACS

EOF
   exit 1
fi


force_update=0
if [ $1 = "-f" ]
then
  echo Forcing update of all files
  force_update=1
  shift
fi

reverse=0
if [ $1 = "-R" ]
then
  echo Reverse: updating git tree from downstream tree
  reverse=1
  shift
fi

# infer source path from name of script
source=$(dirname "$0")

# check general validity of target path
target="$1"
if [ ! -d "${target}" ]
then
    echo ERROR: Target directory ${target} does not exist
    exit 2
fi

# undocumented option to only compare trees
checkonly=0
[ "$2" = "--diff" ] && checkonly=1
[ $force_update = 1 ] && checkonly=0

# try to determine what code resides inside the target dir
code=unkown
if [ -f "${target}/src/lammps.h" ]
then
  code=LAMMPS
elif [ -f "${target}/src/NamdTypes.h" ]
then
  code=NAMD
elif [ -f "${target}/src/VMDApp.h" ]
then
  code=VMD
elif [ -f "${target}/src/gromacs/commandline.h" ]
then
  code=GROMACS
else
  # handle the case if the user points to ${target}/src
  target=$(dirname "${target}")
  if [ -f "${target}/src/lammps.h" ]
  then
    code=LAMMPS
  elif [ -f "${target}/src/NamdTypes.h" ]
  then
    code=NAMD
  elif [ -f "${target}/src/VMDApp.h" ]
  then
    code=VMD
  elif [ -f "${target}/src/gromacs/commandline.h" ]
  then
    code=GROMACS
  else
    echo ERROR: Cannot detect a supported code in the target directory
    exit 3
  fi
fi

echo Detected ${code} source tree in ${target}
echo -n Updating

# conditional file copy
condcopy () {
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
      cmp -s "$a" "$b" || cp "$a" "$b"
      echo -n '.'
    fi
  fi
}


# check files related to, but not part of the colvars module
checkfile () {
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
    echo "Differences found between ${a} and ${b} -- Check $(basename ${a}).diff"
    if [ $force_update = 1 ]
    then
      echo "Overwriting ${b}, as requested by the -f flag."
      cp "$a" "$b"
    fi
  else
    rm -f $(basename ${a}).diff
  fi
}

# update LAMMPS tree
if [ ${code} = LAMMPS ]
then

  # update code-independent headers
  for src in ${source}/src/*.h
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/lib/colvars/${tgt}"
  done
  # update code-independent sources
  for src in ${source}/src/*.cpp
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/lib/colvars/${tgt}"
  done

  # update LAMMPS interface files (library part)
  for src in ${source}/lammps/lib/colvars/Makefile.* ${source}/lammps/lib/colvars/README
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/lib/colvars/${tgt}"
  done
  # update LAMMPS interface files (package part)
  for src in ${source}/lammps/src/USER-COLVARS/*.cpp  ${source}/lammps/src/USER-COLVARS/*.h \
    ${source}/lammps/src/USER-COLVARS/Install.sh ${source}/lammps/src/USER-COLVARS/README
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/src/USER-COLVARS/${tgt}"
  done

  # update LAMMPS documentation
  for src in ${source}/lammps/doc/*.txt
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/doc/${tgt}"
  done

  for src in ${source}/doc/colvars-refman-lammps.pdf
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/doc/PDF/${tgt}"
  done

  echo ' done.'
  exit 0
fi

# update NAMD tree
if [ ${code} = NAMD ]
then

  # update code-independent headers
  for src in ${source}/src/*.h
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/src/${tgt}"
  done
  # update code-independent sources
  for src in ${source}/src/*.cpp
  do \
    tgt=$(basename ${src%.cpp})
    condcopy "${src}" "${target}/src/${tgt}.C"
  done

  # update NAMD interface files
  for src in ${source}/namd/src/colvarproxy_namd.h ${source}/namd/src/colvarproxy_namd.C
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/src/${tgt}"
  done

  condcopy "${source}/doc/colvars-refman.bib" "${target}/ug/ug_colvars.bib"
  condcopy "${source}/doc/colvars-refman-main.tex" "${target}/ug/ug_colvars.tex"
  condcopy "${source}/namd/ug/ug_colvars_macros.tex" "${target}/ug/ug_colvars_macros.tex"
  condcopy "${source}/doc/colvars_diagram.pdf" "${target}/ug/figures/colvars_diagram.pdf"
  condcopy "${source}/doc/colvars_diagram.eps" "${target}/ug/figures/colvars_diagram.eps"

  echo ' done.'

  for src in ${source}/namd/src/* 
  do \
    tgt=$(basename ${src})
    checkfile "${src}" "${target}/src/${tgt}"
  done
  for src in ${source}/namd/Make*
  do 
    tgt=$(basename ${src})
    checkfile "${src}" "${target}/${tgt}"
  done

  exit 0
fi


# update VMD tree
if [ ${code} = VMD ]
then

  # update code-independent headers
  for src in ${source}/src/*.h
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/src/${tgt}"
  done
  # update code-independent sources
  for src in ${source}/src/*.cpp
  do \
    tgt=$(basename ${src%.cpp})
    condcopy "${src}" "${target}/src/${tgt}.C"
  done

  condcopy "${source}/doc/colvars-refman.bib" "${target}/doc/ug_colvars.bib"
  condcopy "${source}/doc/colvars-refman-main.tex" "${target}/doc/ug_colvars.tex"
  condcopy "${source}/vmd/doc/ug_colvars_macros.tex" "${target}/doc/ug_colvars_macros.tex"
  condcopy "${source}/doc/colvars_diagram.pdf" "${target}/doc/pictures/colvars_diagram.pdf"
  condcopy "${source}/doc/colvars_diagram.eps" "${target}/doc/pictures/colvars_diagram.eps"

  # update VMD interface files
  for src in ${source}/vmd/src/colvarproxy_vmd.h ${source}/vmd/src/colvarproxy_vmd.C  
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/src/${tgt}"
  done

  echo ' done.'

  for src in ${source}/vmd/src/* 
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


# update GROMACS tree
if [ ${code} = GROMACS ]
then
  # Copy the colvars source code into gromacs/pulling
  for src in colvaratoms.cpp colvarbias_abf.cpp colvarbias_alb.cpp colvarbias.cpp colvarbias_meta.cpp colvarbias_restraint.cpp colvarcomp_angles.cpp colvarcomp_coordnums.cpp colvarcomp.cpp colvarcomp_distances.cpp colvarcomp_protein.cpp colvarcomp_rotations.cpp colvar.cpp colvargrid.cpp colvarmodule.cpp colvarparse.cpp colvarscript.cpp colvartypes.cpp colvarvalue.cpp colvaratoms.h colvarbias_abf.h colvarbias_alb.h colvarbias.h colvarbias_meta.h colvarbias_restraint.h colvarcomp.h colvargrid.h colvar.h colvarmodule.h colvarparse.h colvarproxy.h colvarscript.h colvartypes.h colvarvalue.h
  do \
    condcopy "src/${src}" "${target}/src/gromacs/pulling/${src}"
  done

  # Copy the GROMACS interface files into gromacs/pulling
  srcDir=${source}/gromacs/src
  for src in colvarproxy_gromacs.cpp colvarproxy_gromacs.h colvars_potential.h
  do \
      condcopy "$srcDir/${src}" "${target}/src/gromacs/pulling/${src}"
  done

  # Find the GROMACS sim_util file, which has changed from
  # sim_util.c to sim_util.cpp between versions.
  if [ -f ${target}/src/gromacs/mdlib/sim_util.cpp ]
  then
      sim_util=${target}/src/gromacs/mdlib/sim_util.cpp
  elif [ -f ${target}/src/gromacs/mdlib/sim_util.c ]
  then
      sim_util=${target}/src/gromacs/mdlib/sim_util.c
  else
      echo "ERROR: Cannot find sim_util.c or sim_util.cpp in the GROMACS source"
      exit 4
  fi
  
  if [ `grep -c colvars $sim_util` -gt 0 ]
  then
      echo "$sim_util appears to already have Colvars modifications. Not modifying."
  else
      # Backup sim_util.
      cp $sim_util ${sim_util}.orig
      
      # Insert necessary pieces of code into the GROMACS sim_util.c or sim_util.cpp.
      awk -f ${source}/gromacs/gromacs-insert.awk ${sim_util}.orig > $sim_util
  fi

  echo ' done.'
  exit 0
fi
