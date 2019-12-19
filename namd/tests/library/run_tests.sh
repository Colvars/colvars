#!/bin/bash
# -*- sh-basic-offset: 2; sh-indentation: 2; -*-

# Run automated tests for NAMD/colvars
# each test is defined by a directory with NAMD input test.namd
# and output files (text only) to be matched in the ExpectedResults/ subdir
# Returns 1 if any test failed, otherwise 0.

# binary to be tested is specified as command-line argument (defaults to namd2)

gen_ref_output=''

export TMPDIR=${TMPDIR:-/tmp}

DIRLIST=''
BINARY=namd2
while [ $# -ge 1 ]; do
  if { echo $1 | grep -q namd2 ; }; then
    echo "Using NAMD executable from $1"
    BINARY=$1
  elif [ "x$1" = 'x-g' ]; then
    gen_ref_output='yes'
  elif [ "x$1" = 'x-h' ]; then
    echo "Usage: ./run_tests.sh [-h] [-g] [path_to_namd2] [testdir1 [testdir2 ...]]"  >& 2
    echo "    The -g option (re)generates reference outputs in the given directories" >& 2
    echo "    If no executable is given, \"namd2\" is used" >& 2
    echo "    If no directories are given, all matches of [0-9][0-9][0-9]_* are used" >& 2
    echo "    This script relies on the executable spiff to be available, and will try to " >& 2
    echo "    download and build it into $TMPDIR if needed." >& 2
    exit 0
  else
    DIRLIST=`echo ${DIRLIST} $1`
  fi
  shift
done

TOPDIR=$(git rev-parse --show-toplevel)
if [ ! -d ${TOPDIR} ] ; then
  echo "Error: cannot identify top project directory." >& 2
  exit 1
fi

SPIFF=$(${TOPDIR}/devel-tools/get_spiff)
if [ $? != 0 ] ; then
    echo "Error: spiff is not available and could not be downloaded/built." >& 2
    exit 1
else
    echo "Using spiff executable from $SPIFF"
    hash -p ${SPIFF} spiff
fi

if ! { echo ${DIRLIST} | grep -q 0 ; } then
  DIRLIST=`eval ls -d [0-9][0-9][0-9]_*`
fi

NUM_THREADS=1
NUM_CPUS=$(nproc)
if [ ${NUM_THREADS} -gt ${NUM_CPUS} ] ; then
  NUM_THREADS=${NUM_CPUS}
fi

TPUT_RED='true'
TPUT_GREEN='true'
TPUT_BLUE='true'
TPUT_CLEAR='true'
if hash tput >& /dev/null ; then
  TPUT_RED='tput setaf 1'
  TPUT_GREEN='tput setaf 2'
  TPUT_BLUE='tput setaf 4'
  TPUT_CLEAR='tput sgr 0'
fi

BASEDIR=$PWD
ALL_SUCCESS=1

# Precision requested to pass (negative powers of ten)
DIFF_PREC=6
# Minimum precision to be tested
MIN_PREC=1

cleanup_files() {
  for script in test*.namd testres*.namd ; do
    if test -L ${script} ; then
      rm -f ${script}
    fi
    for f in ${script%.namd}.*diff; do if [ ! -s $f ]; then rm -f $f; fi; done # remove empty diffs only
    rm -f ${script%.namd}.*{BAK,old,backup}
    for f in ${script%.namd}.*{state,state.stripped,out,traj,coor,vel,xsc,dcd,pmf,hills,grad,force,count,histogram?.dat,hist.dat,corrfunc.dat,histogram?.dx,count.dx,pmf.dx,output.dat}
    do
      if [ ! -f "$f.diff" ]; then rm -f $f; fi # keep files that have a non-empty diff
    done
    rm -f *.out *.out.diff # Delete output files regardless
    rm -f metadynamics1.*.files.txt metadynamics1.*.files.txt.BAK replicas.registry.txt
  done
  tclsh ../Common/delete_tmp_files.tcl
}


for dir in ${DIRLIST} ; do

  if [ -f ${dir}/disabled ] ; then
    continue
  fi

  echo -ne "Entering $(${TPUT_BLUE})${dir}$(${TPUT_CLEAR}) ..."
  cd $dir

  if [ ! -d AutoDiff ] ; then
    echo ""
    echo "  Creating directory AutoDiff, use -g to fill it."
    mkdir AutoDiff
    cd $BASEDIR
    continue
  else

    if [ "x${gen_ref_output}" != 'xyes' ]; then

      if ! { ls AutoDiff/ | grep -q test ; } then
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
    ln -fs ${SCRIPTS} ./
  fi

  # run simulation(s)
  for script in ${SCRIPTS} ; do

    basename=`basename ${script}`
    basename=${basename%.namd}

    # Try running the test (use a subshell to avoid cluttering stdout)
    # Use --source to avoid letting NAMD change its working directory
    # Use multiple threads to test SMP code (TODO: move SMP tests to interface?)
    if ! ( $BINARY +p ${NUM_THREADS} \
                   --source $script > ${basename}.out || \
             false ) > /dev/null 2>&1 ; then
      # This test may be using syntax that changed between versions
      if [ -f ${script%.namd}.legacy.namd ] ; then
        # Try a legacy input
        ( $BINARY +p ${NUM_THREADS} \
                  --source ${script%.namd}.legacy.namd > ${basename}.out || \
            false ) > /dev/null 2>&1
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

    if [ -f ${basename}.colvars.state ] ; then
      # Filter out the version number from the state files to allow comparisons
      grep -sv '^  version ' ${basename}.colvars.state | \
        grep -sv '^  units ' \
        > ${TMPDIR}/${basename}.colvars.state.stripped
      mv -f ${TMPDIR}/${basename}.colvars.state.stripped ${basename}.colvars.state.stripped
    fi

    # If this test is used to generate the reference output files, copy them
    if [ "x${gen_ref_output}" = 'xyes' ]; then
      grep 'NAMD' ${basename}.out | head -n 1 > namd-version.txt
      grep 'Initializing the collective variables module, version' ${basename}.out | head -n 1 >> namd-version.txt
      grep 'Using NAMD interface, version' ${basename}.out | head -n 1 >> namd-version.txt
      cp ${basename}.colvars.state.stripped AutoDiff/
      cp ${basename}.colvars.traj           AutoDiff/
      cp ${basename}.colvars.out            AutoDiff/
      if [ -f ${basename}.histogram1.dat ] ; then
        cp -f ${basename}.histogram1.dat AutoDiff/
      fi
      if [ -f ${basename}.pmf ] ; then
        cp -f ${basename}.pmf AutoDiff/
      fi
    fi

    # Old versions did not accurately update the prefix
    if [ -f .histogram1.dat ] ; then
      mv .histogram1.dat ${basename}.histogram1.dat
    fi

  done

  # now check results
  SUCCESS=1
  for f in AutoDiff/*
  do
    base=`basename $f`
    if [ ! -f $base ] ; then
      echo -e "\n*** File $(${TPUT_RED})$base$(${TPUT_CLEAR}) is missing. ***"
      SUCCESS=0
      ALL_SUCCESS=0
      break
    fi

    if [ "${base}" != "${base%.traj}" ] ; then
      # System force is now total force
      sed 's/fs_/ft_/g' < ${base} > ${TMPDIR}/${base}
      mv -f ${TMPDIR}/${base} ${base}
    fi
    if [ ${base} != ${base%.out} ] ; then
      # Lots of text confuse spiff
      diff $f $base > "$base.diff"
      RETVAL=$?
    else
      ${SPIFF} -r 1e-${DIFF_PREC} $f $base > "$base.diff"
      RETVAL=$?
    fi
    if [ $RETVAL -ne 0 ]
    then
      if [ ${base} != ${base%.out} ]
      then
        echo -n "(warning: differences in log file $base) "
      else
        echo -e "\n*** Failure for file $(${TPUT_RED})$base$(${TPUT_CLEAR}): see `pwd`/$base.diff "
        SUCCESS=0
        ALL_SUCCESS=0
        LOW_PREC=${DIFF_PREC}
        RETVAL=1
        while [ $RETVAL -ne 0 ] && [ $LOW_PREC -gt $MIN_PREC ]
        do
          LOW_PREC=$((${LOW_PREC} - 1))
          spiff -r 1e-${LOW_PREC} $f $base > /dev/null
          RETVAL=$?
        done
        if [ $RETVAL -eq 0 ]
        then
          echo " --> Passes at reduced precision 1e-${LOW_PREC}"
        else
          echo " --> Fails at minimum tested precision 1e-${LOW_PREC}"
        fi
      fi
    fi
  done

  if [ $SUCCESS -eq 1 ]
  then
    if [ "x${gen_ref_output}" == 'xyes' ]; then
      echo "Reference files copied successfully."
    else
      echo "$(${TPUT_GREEN})Success!$(${TPUT_CLEAR})"
    fi
    cleanup_files
  fi

  # TODO: at this point, we may use the diff file to update the reference tests for harmless changes
  # (e.g. keyword echos). Before then, figure out a way to strip the formatting characters produced by spiff.

  cd $BASEDIR
done


if [ $ALL_SUCCESS -eq 1 ]
then
  echo "$(${TPUT_GREEN})All tests succeeded.$(${TPUT_CLEAR})"
  exit 0
else
  echo "$(${TPUT_RED})There were failed tests.$(${TPUT_CLEAR})"
  exit 1
fi
