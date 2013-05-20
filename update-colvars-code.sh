#!/bin/sh
# Script to update a LAMMPS or NAMD source tree with the latest colvars code.

if [ $# -lt 1 ]
then
    cat <<EOF

 usage: sh $0 <target source tree>

   "target source tree" = root directory of the MD code sources
   supported MD codes: LAMMPS, NAMD

EOF
   exit 1
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

# try to determine what code resides inside the target dir
code=unkown
if [ -f "${target}/src/lammps.h" ]
then
  code=LAMMPS
elif [ -f "${target}/src/NamdTypes.h" ]
then
  code=NAMD
else
  # handle the case if the user points to ${target}/src
  target=$(dirname "${target}")
  if [ -f "${target}/src/lammps.h" ]
  then
    code=LAMMPS
  elif [ -f "${target}/src/NamdTypes.h" ]
  then
    code=NAMD
  else
    echo ERROR: Cannot detect a supported code in the target directory
    exit 3
  fi
fi

echo Detected ${code} source tree in ${target}
echo Beginning update...

# conditional file copy
condcopy () {
  if [ -d $(dirname "$2") ]
  then
    cmp -s "$1" "$2" || cp "$1" "$2"
  fi
}

# update LAMMPS tree
if [ ${code} = LAMMPS ]
then
  # update colvars library headers
  for src in ${source}/src/*.h
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/lib/colvars/${tgt}"
  done
  # update colvars library code
  for src in ${source}/src/*.C
  do \
    tgt=$(basename ${src%.C})
    condcopy "${src}" "${target}/lib/colvars/${tgt}.cpp"
  done
  # update LAMMPS infrastructure for colvars library
  for src in ${source}/lammps/lib/colvars/*
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/lib/colvars/${tgt}"
  done
  # update LAMMPS interface code for colvars
  for src in ${source}/lammps/src/USER-COLVARS/*
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/src/USER-COLVARS/${tgt}"
  done
  echo Update complete.
  exit 0
fi

# update NAMD tree
if [ ${code} = NAMD ]
then
  # update colvars library
  for src in ${source}/src/*.h ${source}/src/*.C
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/src/${tgt}"
  done
  # update NAMD interface files
  for src in ${source}/namd/src/*
  do \
    tgt=$(basename ${src})
    condcopy "${src}" "${target}/src/${tgt}"
  done
  echo Update complete.
  exit 0
fi
