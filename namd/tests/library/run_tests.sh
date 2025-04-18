#!/bin/bash
# -*- sh-basic-offset: 2; sh-indentation: 2; -*-

# Run automated tests for NAMD/colvars
# each test is defined by a directory with NAMD input test.namd
# and output files (text only) to be matched in the ExpectedResults/ subdir
# Returns 1 if any test failed, otherwise 0.

# binary to be tested is specified as command-line argument (defaults to namd3)

gen_ref_output=''

export TMPDIR=${TMPDIR:-/tmp}

DIRLIST=''
BINARY=namd3
CUDASOA=0

while [ $# -ge 1 ]; do
  if { echo $1 | grep -q namd ; }; then
    echo "Using NAMD executable from $1"
    BINARY=$1
  elif [ "x$1" = 'x-g' ]; then
    gen_ref_output='yes'
    echo "Generating reference output"
  elif [ "x$1" = 'x-h' ]; then
    echo "Usage: ./run_tests.sh [-h] [-g] [-cudasoa] [path_to_namd] [testdir1 [testdir2 ...]]"  >& 2
    echo "    The -g option (re)generates reference outputs in the given directories" >& 2
    echo "    The -cudasoa option enables the GPU-resident NAMD3 code path (CUDASOAIntegrate)" >& 2
    echo "    If no executable is given, \"namd3\" is used" >& 2
    echo "    If no directories are given, all matches of [0-9][0-9][0-9]_* are used" >& 2
    echo "    This script relies on the executable spiff to be available, and will try to " >& 2
    echo "    download and build it into $TMPDIR if needed." >& 2
    exit 0
  elif [ "x$1" = 'x-cudasoa' ]; then
    CUDASOA=1
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


NUM_THREADS=4
NUM_TASKS=1

CHARM_ARCH=$(${BINARY} 2>&1 | grep 'Info: Based on Charm++/Converse' | cut -d' ' -f 7 || true)
if [ "x${CHARM_ARCH}" == "xmpi-linux-x86_64" ] && source ${TOPDIR}/devel-tools/load-openmpi.sh ; then
  NUM_TASKS=${NUM_THREADS}
  NUM_THREADS=1
  CMD="mpirun -n ${NUM_TASKS} -oversubscribe $BINARY"
else
  CMD="$BINARY +p ${NUM_THREADS}"
fi

echo "Running NAMD as: $CMD"

TPUT_RED='true'
TPUT_GREEN='true'
TPUT_YELLOW='true'
TPUT_BLUE='true'
TPUT_CLEAR='true'
if hash tput >& /dev/null && [ -z "${GITHUB_ACTION}" ] ; then
  TPUT_RED='tput setaf 1'
  TPUT_GREEN='tput setaf 2'
  TPUT_YELLOW='tput setaf 3'
  TPUT_BLUE='tput setaf 4'
  TPUT_CLEAR='tput sgr 0'
fi

BASEDIR=$PWD
ALL_SUCCESS=1

# Precision requested to pass (negative powers of ten)
DIFF_PREC=6
DIFF_ABS_PREC=12
# Minimum precision to be tested
MIN_PREC=1

cleanup_files() {
  for script in test*.namd testres*.namd ; do
    if test -L ${script} ; then
      rm -f ${script}
    fi
    for f in ${script%.namd}.*diff; do if [ ! -s $f ]; then rm -f $f; fi; done # remove empty diffs only
    rm -f ${script%.namd}.*{BAK,old,backup}
    for f in ${script%.namd}.*{state,state.stripped,out,traj,coor,vel,xsc,dcd,pmf,hills,grad,force,count,histogram?.dat,hist.dat,corrfunc.dat,histogram?.dx,count.dx,pmf.dx,output.dat,ti,kernels.dat}
    do
      if [ ! -f "$f.diff" ]; then rm -f $f; fi # keep files that have a non-empty diff
    done
    rm -f *.out *.out.diff # Delete output files regardless
    rm -f metadynamics1.*.files.txt metadynamics1.*.files.txt.BAK replicas.registry.txt
  done
  tclsh ../Common/delete_tmp_files.tcl
}

declare -a failed_tests
declare -a failed_tests_low_prec
declare -a failed_tests_abs_prec

TORCH_LINKED=false
if { ldd $(which $BINARY) | grep -q libtorch[_a-zA-Z]*.so ; } then TORCH_LINKED=true ; fi

for dir in ${DIRLIST} ; do

  if [ -f ${dir}/disabled ] ; then
    continue
  fi

  if [ -f ${dir}/skip_test.sh ]; then
    bash ${dir}/skip_test.sh
    if [ $? -eq 0 ]; then
      echo "Directory ${dir} skipped."
      continue
    fi
  fi

  if echo ${dir} | grep -q torchann && [ ${TORCH_LINKED} != "true" ] ; then 
    echo "Directory ${dir} skipped."
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
  fi

  # run simulation(s)
  for script in ${SCRIPTS} ; do

    if [ ! -f `basename ${script}` ] ; then
      ln -fs ${script} ./
    fi

    script=`basename ${script}`
    basename=${script%.namd}

    # If we are doing binary restarts, make an exception for tests that don't support it
    if [ -n "${COLVARS_BINARY_RESTART}" ] ; then
      if [ ${dir} == 004_10ala_moving_restart ] ; then
        export COLVARS_BINARY_RESTART=0
      else
        export COLVARS_BINARY_RESTART=1
      fi
    fi

    if [ -f no_smp ]; then
      # Force the test to run serially
      NAMD_CUDASOA=$CUDASOA $BINARY +p1 $script > ${basename}.out
    else
      # Run the test (use a subshell to avoid cluttering stdout)
      NAMD_CUDASOA=$CUDASOA $CMD $script > ${basename}.out
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
      echo -n " (Copying reference files for ${script}) "
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
      if [ -f ${basename}.grad ] ; then
        cp -f ${basename}.grad AutoDiff/
        cp -f ${basename}.count AutoDiff/
      fi

      # Update any additional files with current versions
      for file in AutoDiff/${basename}*; do
        source=`basename ${file}`
        if [ -f ${source} ] ; then
          cp -uf ${source} AutoDiff/
        fi
      done
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
    if [ "${base%.state.stripped}" != "${base}" ] && [ -n "${COLVARS_BINARY_RESTART}" ] ; then
      # Do not try comparing binary state files, they will never match anyway
      continue
    fi
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
      if [ ${base} == ${base%.out} ] # Ignore differences in stdout log
      then
        echo -e "\n*** Failure for file $(${TPUT_RED})$base$(${TPUT_CLEAR}): see `pwd`/$base.diff "
        SUCCESS=0
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
          ALL_SUCCESS=0
          failed_tests_low_prec+=($dir)
          echo " --> Passes at reduced precision 1e-${LOW_PREC}"
        else
	  # Test absolute error
	  spiff -a 1e-${DIFF_ABS_PREC} $f $base > /dev/null
	  RETVAL=$?
	  if [ $RETVAL -eq 0 ]
	  then
	    failed_tests_abs_prec+=($dir)
	    echo " --> Passes at absolute precision 1e-${DIFF_ABS_PREC}"
	  else
            ALL_SUCCESS=0
            failed_tests+=($dir)
            echo " --> Fails at minimum tested precision 1e-${LOW_PREC}"
          fi
        fi
      fi
    fi
  done

  if [ $SUCCESS -eq 1 ]
  then
    if [ "x${gen_ref_output}" == 'xyes' ]; then
      echo "Reference files copied successfully."
    else
      echo " $(${TPUT_GREEN})Success!$(${TPUT_CLEAR})"
    fi
    cleanup_files
  fi

  # TODO: at this point, we may use the diff file to update the reference tests for harmless changes
  # (e.g. keyword echos). Before then, figure out a way to strip the formatting characters produced by spiff.

  cd $BASEDIR
done

if [ ${#failed_tests_abs_prec[@]} -gt 0 ]; then
  echo "$(${TPUT_YELLOW})The following tests are failed at relative precision, but passed at absolute precision:$(${TPUT_CLEAR})"
  printf "%s\n" "${failed_tests_abs_prec[@]}" | sort -u
fi

if [ $ALL_SUCCESS -eq 1 ]
then
  echo "$(${TPUT_GREEN})All tests succeeded.$(${TPUT_CLEAR})"
  exit 0
else
  echo "$(${TPUT_RED})There were failed tests.$(${TPUT_CLEAR})"
  if [ ${#failed_tests[@]} -gt 0 ]; then
    echo "The following tests are failed:"
    printf "%s\n" "${failed_tests[@]}" | sort -u
  fi
  if [ ${#failed_tests_low_prec[@]} -gt 0 ]; then
    echo "The following tests are failed, but passed at low precisions:"
    printf "%s\n" "${failed_tests_low_prec[@]}" | sort -u
  fi
  exit 1
fi
