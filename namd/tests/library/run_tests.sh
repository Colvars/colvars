#!/bin/bash
# -*- sh-basic-offset: 2; sh-indentation: 2; -*-

# Run automated tests for NAMD/colvars
# each test is defined by a directory with NAMD input test.namd
# and output files (text only) to be matched in the ExpectedResults/ subdir
# Returns 1 if any test failed, otherwise 0.

# binary to be tested is specified as command-line argument (defaults to namd2)

gen_ref_output=''

DIRLIST=''
BINARY=namd2
while [ $# -ge 1 ]; do
  if { echo $1 | grep -q namd2 ; }; then
    BINARY=$1
  elif [ "x$1" = 'x-g' ]; then
    gen_ref_output='yes'
  else
    DIRLIST=`echo ${DIRLIST} $1`
  fi
  shift
done
if ! { echo ${DIRLIST} | grep -q 0 ; } then
  DIRLIST=`eval ls -d [0-9][0-9][0-9]_*`
fi

DIFF=spiff
BASEDIR=$PWD
ALL_SUCCESS=1


cleanup_files() {
  for script in test*.namd testres*.namd ; do
    for f in ${script%.namd}.*diff; do if [ ! -s $f ]; then rm -f $f; fi; done # remove empty diffs only
    rm -f ${script%.namd}.*{BAK,old,backup}
    for f in ${script%.namd}.*{state,out,traj,coor,vel,xsc,pmf,hills,grad,count}
    do
      if [ ! -f "$f.diff" ]; then rm -f $f; fi # keep files that have a non-empty diff
    done
    rm -f metadynamics1.*.files.txt replicas.registry.txt
  done
}


for dir in ${DIRLIST} ; do
  echo -ne "Entering $dir ..."
  cd $dir

  if [ ! -d AutoDiff ] ; then
    echo ""
    echo "  Creating directory AutoDiff, use -g to fill it."
    mkdir AutoDiff
    cd $BASEDIR
    continue
  else

    if [ "x${gen_ref_output}" != 'xyes' ]; then

      if ! { ls AutoDiff/ | grep -q traj ; } then
        echo ""
        echo "  Warning: directory AutoDiff empty!"
        cd $BASEDIR
        continue
      fi

      # first, remove target files from work directory
      for f in AutoDiff/*
      do
        base=`basename $f`
        if [ -f $base ]
        then
          mv $base $base.backup
        fi
      done
    fi
  fi

  cleanup_files

  if ls | grep -q \.namd ; then
    SCRIPTS=`ls -1 *namd | grep -v legacy`
  else
    SCRIPTS="../Common/test.namd ../Common/test.restart.namd"
  fi

  # run simulation(s)
  for script in ${SCRIPTS} ; do

    basename=`basename ${script}`
    basename=${basename%.namd}

    # Try running the test (use a subshell to avoid cluttering stdout)
    # Use --source to avoid letting NAMD change its working directory
    # Use multiple threads to test SMP code (TODO: move SMP tests to interface?)
    if ! ( $BINARY +p 3 --source $script > ${basename}.out || false ) > /dev/null 2>&1 ; then
      # This test may be using syntax that changed between versions
      if [ -f ${script%.namd}.legacy.namd ] ; then
        # Try a legacy input
        ( $BINARY +p 3 --source ${script%.namd}.legacy.namd > ${basename}.out || false ) > /dev/null 2>&1
      fi
    fi

    # Output of Colvars module, minus the version numbers
    grep "^colvars:" ${basename}.out | grep -v 'Initializing the collective variables module' \
      | grep -v 'Using NAMD interface, version' > ${basename}.colvars.out

    # Output of Tcl interpreter for automatic testing of scripts (TODO: move this to interface)
    grep "^TCL:" ${basename}.out | grep -v '^TCL: Suspending until startup complete.' > ${basename}.Tcl.out
    if [ ! -s ${basename}.Tcl.out ]; then
      rm -f ${basename}.Tcl.out
    fi

    # If this test is used to generate the reference output files, copy them
    if [ "x${gen_ref_output}" = 'xyes' ]; then
      cp ${basename}.colvars.state AutoDiff/
      cp ${basename}.colvars.traj  AutoDiff/
      cp ${basename}.colvars.out   AutoDiff/
    fi

  done

  # now check results
  SUCCESS=1
  for f in AutoDiff/*
  do
    base=`basename $f`
    if [ "${base}" != "${base%.state}" ] ; then
      # Filter out the version number from the state files to allow comparisons
      grep -v 'version' ${base} > ${base}.tmp
      mv ${base}.tmp ${base}
    fi
    if [ "${base}" != "${base%.traj}" ] ; then
      # System force is now total force
      sed 's/fs_/ft_/g' < ${base} > ${base}.tmp
      mv ${base}.tmp ${base}
    fi
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
    if [ "x${gen_ref_output}" == 'xyes' ]; then
      echo "Reference files copied successfully."
    else
      echo "Success!"
    fi
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
