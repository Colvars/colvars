#!/bin/bash
# -*- sh-basic-offset: 2; sh-indentation: 2; -*-

# Run automated tests for NAMD/colvars
# each test is defined by a directory with NAMD input test.namd
# and output files (text only) to be matched in the ExpectedResults/ subdir
# Returns 1 if any test failed, otherwise 0.

# binary to be tested is specified as command-line argument (defaults to namd2)



BINARY=namd2
if [ $# -ge 1 ]; then
  if { echo $1 | grep -q namd2 ; }; then
    BINARY=$1
    shift
  fi
fi
DIRLIST=`eval ls -d [0-9][0-9][0-9]_*`
if [ $# -ge 1 ]; then
  DIRLIST=`echo $@`
fi

DIFF=spiff
BASEDIR=$PWD
ALL_SUCCESS=1


cleanup_files() {
  for script in test*.namd testres*.namd ; do
    for f in ${script%.namd}.*diff; do if [ ! -s $f ]; then rm -f $f; fi; done # remove empty diffs only
    rm -f ${script%.namd}.*{BAK,old,backup}
    rm -f ${script%.namd}.*{state,out,traj,coor,vel,xsc,pmf,hills,grad,count}
    rm -f metadynamics1.*.files.txt replicas.registry.txt
  done
}


for dir in ${DIRLIST} ; do
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

  cleanup_files

  # run simulation(s)
  for script in test*.namd ; do
    # use --source to avoid letting NAMD change directory
    # use 5 threads to test SMP code
    $BINARY +p 5 --source $script > ${script%.namd}.out
    # collect output of colvars module, except the version numbers
    grep "^colvars:" ${script%.namd}.out | grep -v 'Initializing the collective variables module' \
      | grep -v 'Using NAMD interface, version' > ${script%.namd}.colvars.out
    # Output of Tcl interpreter for automatic testing of scripts
    grep "^TCL:" ${script%.namd}.out | grep -v '^TCL: Suspending until startup complete.' > ${script%.namd}.Tcl.out
    if [ ! -s ${script%.namd}.Tcl.out ]; then
      rm -f ${script%.namd}.Tcl.out
    fi

    # Filter out the version number from the state files to allow comparisons
    grep -v 'version' ${script%.namd}.colvars.state > ${script%.namd}.colvars.state.tmp
    mv ${script%.namd}.colvars.state.tmp ${script%.namd}.colvars.state

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
        echo -e "\n***  Failure for file $base: see `pwd`/$base.diff ***"
        SUCCESS=0
        ALL_SUCCESS=0
      fi
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

  for dir in [0-9][0-9][0-9]_* ; do
    cd $dir
    cleanup_files
    cd $BASEDIR
  done
  
  exit 0
else
  echo "There were failed tests."
  exit 1
fi

