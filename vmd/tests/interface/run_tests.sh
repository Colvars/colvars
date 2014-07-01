#!/bin/bash

# Run automated tests for VMD/colvars
# each test is defined by a directory with VMD input test.tcl
# and output files (text only) to be matched in the ExpectedResults/ subdir
# Returns 1 if any test failed, otherwise 0.

# binary to be tested is specified as command-line argument (defaults to namd2)

if [ $# -lt 1 ]
then
  BINARY=vmd
else
  BINARY=$1
fi

DIFF=spiff
BASEDIR=$PWD
ALL_SUCCESS=1
TMP=/tmp


cleanup_files () {
for dir in [0-9][0-9][0-9]_* ; do
  for script in test*.vmd ; do
       rm -f ${dir}/${script%.vmd}.*{state,out,traj,coor,vel,xsc,BAK,old,backup,diff,pmf,grad,count}
  done
done
}



for dir in [0-9][0-9][0-9]_*
do
  echo -ne "Entering $dir... "
  cd $dir

  # first, remove target files from work directory
  for f in AutoDiff/*
  do
    base=`basename $f`
    if [ -f $base ]
    then
      mv $base $base.backup
    fi
  done

  # run simulation(s)
  for script in test*.vmd ; do
      $BINARY -dispdev none -e $script > ${script%.vmd}.out
      # collect output of colvars module, except the version number
      # TODO: strip also the echo of newly introduced keywords
      grep "^colvars:" ${script%.vmd}.out | grep -v 'Initializing the collective variables module' > ${script%.vmd}.colvars.out
  done
  
  # now check results
  SUCCESS=1
  for f in AutoDiff/*
  do
    base=`basename $f`
    $DIFF $f $base > "$base.diff"
    RETVAL=$?
    if [ $RETVAL -ne 0 ]
    then
      echo "***  Failure for file $base: see $dir/$base.diff ***"
      SUCCESS=0
      ALL_SUCCESS=0
    fi
  done

  if [ $SUCCESS -eq 1 ]
  then
    echo "Success!"
  fi

  cd $BASEDIR
done


if [ $ALL_SUCCESS -eq 1 ]
then
  echo "All tests succeeded."
  cleanup_files
  exit 0
else
  echo "There were failed tests."
  exit 1
fi

