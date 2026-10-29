
shopt -s nullglob

TMPDIR=/tmp
DIRLIST=''
SPIFF=spiff

cleanup='true'
verbose='false'
gen_ref_output='false'

while [ $# -ge 1 ]; do
  if [ "x$1" = 'x-g' ]; then
    gen_ref_output='true'
    echo "Generating reference output"
  elif [ "x$1" = 'x-v' ]; then
    verbose='true'
  elif [ "x$1" = 'x-k' ]; then
    cleanup='false'
  elif [ "x$1" = 'x-h' ]; then
    echo "Usage: ./run_tests.sh [-h] [-g] [-v] [-k] [path_to_gmx] [testdir1 [testdir2 ...]]"  >& 2
    echo >& 2
    echo "        If no gmx executable is given, \"gmx_d\" is used" >& 2
    echo >& 2
    echo "        If no directories are given, all matches of [0-9][0-9][0-9]_* are used" >& 2
    echo >& 2
    echo "        This script relies on the executable spiff to be available, and will try to " >& 2
    echo "        download and build it into $TMPDIR if needed." >& 2
    echo >& 2
    echo "  -g    (re)generate reference outputs in the given directories" >& 2
    echo >& 2
    echo "  -k    keep the output files after a successful test, otherwise they are deleted." >& 2
    echo "  -v    enable verbose output." >& 2
    echo >& 2
    echo "        The OpenMM platform is chosen by OpenMM, and can be set with the" >& 2
    echo "        environment variable OPENMM_DEFAULT_PLATFORM (e.g. Reference, CPU, CUDA)." >& 2
    exit 0
  else
    DIRLIST=`echo ${DIRLIST} $1`
  fi
  shift
done

# If no test is specified, run all tests in this directory
#
if [ -z "$DIRLIST" ]; then
    DIRLIST=`eval ls -d [0-9][0-9][0-9]_*`
fi


TPUT_RED='true'
TPUT_GREEN='true'
TPUT_BLUE='true'
TPUT_CLEAR='true'
if hash tput >& /dev/null && [ -z "${GITHUB_ACTION}" ]; then
  TPUT_RED='tput setaf 1'
  TPUT_GREEN='tput setaf 2'
  TPUT_BLUE='tput setaf 4'
  TPUT_CLEAR='tput sgr 0'
fi

BASEDIR=$PWD
ALL_SUCCESS=1

# Precision requested to pass (negative powers of ten)
DIFF_PREC=4
# Minimum precision to be tested
MIN_PREC=1

cleanup_files() {
  for script in test*.in ; do
    for f in ${script%.in}.*diff; do if [ ! -s $f ]; then rm -f $f; fi; done # remove empty diffs only
    rm -f ${script%.in}.*{BAK,old,backup}
    for f in ${script%.in}.*{state,state.stripped,lmp.data,out,traj,pmf,hills,grad,count,{histogram?,corrfunc,hist}.dat,.histogram?.dx,.ti.force}
    do
      if [ ! -f "$f.diff" ]; then rm -f $f; fi # keep files that have a non-empty diff
    done
    rm -f metadynamics1.*.files.txt replicas.registry.txt rest.colvars.state* *.log
  done
}


for dir in $DIRLIST; do

  echo -ne "Entering $(${TPUT_BLUE})${dir}$(${TPUT_CLEAR}) ... "
  if ! cd "$dir"; then
    FAILED_TESTS="$FAILED_TESTS $dir "
    ALL_SUCCESS=0
    continue
  fi

  if [ ! -d AutoDiff ] ; then
    echo ""
    echo "  Creating directory AutoDiff"
    mkdir AutoDiff
    # Still run the tests to check for successful completion
  fi

  if $gen_ref_output; then

    if ! { ls AutoDiff/ | grep -q test ; } then
      echo "  Directory AutoDiff is empty, populating..."
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

  if $cleanup; then
    cleanup_files
  fi
  basename="test"

  ln -fs ../Common/index.ndx
  pixi run --manifest-path ../../pixi.toml python ../runner.py > $basename.log 2>&1

  if [ $? -ne 0 ]; then
    echo "Test $dir failed. See $basename.log for details."
    FAILED_TESTS="$FAILED_TESTS $dir "
    SUCCESS=0
    ALL_SUCCESS=0
    cd "$BASEDIR"
    continue
  else
    SUCCESS=1
  fi

  logfile=${basename}.log
  # Printed by runner.py; may be interleaved with the Colvars output
  echo -n "[$(grep -m 1 -o 'OpenMM platform: .*' ${logfile} | sed 's/^OpenMM platform: //')] "
  # Filter out the version numbers to allow comparisons
  grep "^colvars:" ${logfile} \
  | grep -v 'Initializing the collective variables module' \
  | grep -v 'Using OpenMM interface' \
  | grep -v 'OpenMM platform:' > ${basename}.colvars.out

  if [ -s ${logfile%.log}.colvars.state ] ; then
    grep -sv 'version' ${logfile%.log}.colvars.state \
      > ${TMPDIR}/${logfile%.log}.colvars.state.stripped && \
    mv -f ${TMPDIR}/${logfile%.log}.colvars.state.stripped ${logfile%.log}.colvars.state.stripped
  fi

  # If this test is used to generate the reference output files, copy them
  if $gen_ref_output; then
    grep 'Initializing the collective variables module, version' ${logfile} | head -n 1 >> openmm-colvars-version.txt
    [ -f ${basename}.colvars.state.stripped ] && cp ${basename}.colvars.state.stripped AutoDiff/
    [ -f ${basename}.colvars.traj ] && cp ${basename}.colvars.traj           AutoDiff/
    [ -f ${basename}.colvars.out ] && cp ${basename}.colvars.out            AutoDiff/
    if [ -f ${basename}.histogram1.dat ] ; then
      cp -f ${basename}.histogram1.dat AutoDiff/
    fi
    if [ -f ${basename}.pmf ] ; then
      cp -f ${basename}.pmf AutoDiff/
    fi
  fi

  # now check results
  SUCCESS=1
  for f in AutoDiff/*
  do
    base=`basename $f`
    if $verbose ; then
      echo -ne "  Comparing $(${TPUT_BLUE})$base$(${TPUT_CLEAR}) ..."
    fi

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
        echo -e "\n*** Failure for file $(${TPUT_RED})$base$(${TPUT_CLEAR}) (return code = $RETVAL): see also `pwd`/$base.diff "
        SUCCESS=0
        ALL_SUCCESS=0
        LOW_PREC=${DIFF_PREC}
        RETVAL=1
        while [ $RETVAL -ne 0 ] && [ $LOW_PREC -gt $MIN_PREC ]
        do
          LOW_PREC=$((${LOW_PREC} - 1))
          ${SPIFF} -r 1e-${LOW_PREC} $f $base > /dev/null
          RETVAL=$?
        done
        if [ $RETVAL -eq 0 ]
        then
          echo " --> Passes at reduced precision 1e-${LOW_PREC}"
        else
          echo " --> Fails at minimum tested precision 1e-${LOW_PREC}"
        fi
      elif $verbose ; then
        echo -e "Skipping comparison for file $(${TPUT_BLUE})$base$(${TPUT_CLEAR}) (stdout log)"
      fi
    elif $verbose ; then
      echo -e "Successful comparison for file $(${TPUT_GREEN})$base$(${TPUT_CLEAR})"
    fi
  done

  if [ $SUCCESS -eq 1 ]
  then
    if $gen_ref_output; then
      echo "Reference files copied successfully."
    else
      echo " $(${TPUT_GREEN})Success!$(${TPUT_CLEAR})"
    fi
    if $cleanup; then
      cleanup_files
    fi
  else
    FAILED_TESTS="$FAILED_TESTS $dir "
  fi
  cd "$BASEDIR"
done

if [ -n "$FAILED_TESTS" ]; then
    echo "The following tests failed: $FAILED_TESTS"
    exit 1
fi
