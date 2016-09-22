#!/bin/bash

# Run automated tests for VMD/colvars
# each test is defined by a directory with VMD input test.tcl
# and output files (text only) to be matched in the ExpectedResults/ subdir
# Returns 1 if any test failed, otherwise 0.

# binary to be tested is specified as command-line argument (defaults to vmd)

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
  for script in test*.vmd ; do
    for f in ${script%.vmd}.*diff; do if [ ! -s $f ]; then rm -f $f; fi; done # remove empty diffs only
    for f in ${script%.vmd}.*{state,out,err,traj,coor,vel,xsc,BAK,old,backup,pmf,grad,count}
    do
      if [ ! -f "$f.diff" ]; then rm -f $f; fi # keep files that have a non-empty diff 
    done
  done
}



for dir in [0-9][0-9][0-9]_*
do
  echo "Entering $dir... "
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
      grep "^colvars:" ${script%.vmd}.out | grep -v 'Initializing the collective variables module' | grep -v 'Using VMD interface' > ${script%.vmd}.colvars.out
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
      if [ ${base##*\.} = 'out' ]
      then
        echo -n "(warning: differences in log file $base) "
      else
        echo -e "\n***  Failure for file $base: see $dir/$base.diff ***"
        SUCCESS=0
        ALL_SUCCESS=0
      fi
    fi
  done

  if [ $SUCCESS -eq 1 ]
  then
    echo "Success!"
    cleanup_files
  fi

  cd $BASEDIR
done


if [ $ALL_SUCCESS -eq 1 ]
then
  echo "All tests succeeded."
  exit 0
else
  echo "There were failed tests."
  exit 1
fi

