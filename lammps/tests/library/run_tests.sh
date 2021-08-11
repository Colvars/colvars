#!/bin/bash
# -*- sh-basic-offset: 2; sh-indentation: 2; -*-

# Run automated tests for LAMMPS/Colvars
gen_ref_output=''

export TMPDIR=${TMPDIR:-/tmp}

DIRLIST=''
BINARY=lmp_cpu
while [ $# -ge 1 ]; do
  if [ -x $1 ] && { echo $1 | grep -q lmp ; }; then
    echo "Using LAMMPS executable from $1"
    BINARY=$1
  elif [ "x$1" = 'x-g' ]; then
    gen_ref_output='yes'
  elif [ "x$1" = 'x-h' ]; then
    echo "Usage: ./run_tests.sh [-h] [-g] [path_to_lammps] [testdir1 [testdir2 ...]]"  >& 2
    echo "    The -g option (re)generates reference outputs in the given directories" >& 2
    echo "    If no executable is given, \"lmp_cpu\" is used" >& 2
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
  for script in test*.lmp.in testres*.lmp.in ; do
    for f in ${script%.lmp.in}.*diff; do if [ ! -s $f ]; then rm -f $f; fi; done # remove empty diffs only
    rm -f ${script%.lmp.in}.*{BAK,old,backup}
    for f in ${script%.lmp.in}.*{state,state.stripped,lmp.data,out,traj,pmf,hills,grad,count,histogram?.dat,histogram?.dx}
    do
      if [ ! -f "$f.diff" ]; then rm -f $f; fi # keep files that have a non-empty diff
    done
    rm -f metadynamics1.*.files.txt replicas.registry.txt rest.colvars.state* log.cite log.lammps
  done
}


for dir in ${DIRLIST} ; do
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

  if ls | grep -q \.lmp.in ; then
    SCRIPTS=`ls -1 *lmp.in`
  else
    SCRIPTS="../common/test.lmp.in ../common/test.restart.lmp.in"
  fi

  # run simulation(s)
  for script in ${SCRIPTS} ; do

    basename=`basename ${script}`
    basename=${basename%.lmp.in}

    colvars_config=test.in
    if [ ${script} = ../common/test.restart.lmp.in ] ; then
      if [ -f test.restart.in ] ; then
        colvars_config=test.restart.in
      fi
    fi

    $BINARY -in $script -var colvars_config ${colvars_config} \
        -log ${basename}.out > /dev/null

    # Output of Colvars module, minus the version numbers
    grep "^colvars:" ${basename}.out | grep -v 'Initializing the collective variables module' \
      | grep -v 'Using NAMD interface, version' > ${basename}.colvars.out

    # # Output of Tcl interpreter for automatic testing of scripts (TODO: move this to interface)
    # grep "^TCL:" ${basename}.out | grep -v '^TCL: Suspending until startup complete.' > ${basename}.Tcl.out
    # if [ ! -s ${basename}.Tcl.out ]; then
    #   rm -f ${basename}.Tcl.out
    # fi

    if [ -f ${basename}.colvars.state ] ; then
      # Filter out the version number from the state files to allow comparisons
      grep -sv '^  version' ${basename}.colvars.state | \
        grep -sv '^  units' \
        > ${TMPDIR}/${basename}.colvars.state.stripped
      mv -f ${TMPDIR}/${basename}.colvars.state.stripped ${basename}.colvars.state.stripped
    fi

    # If this test is used to generate the reference output files, copy them
    if [ "x${gen_ref_output}" = 'xyes' ]; then
      head -n 1 ${basename}.out > lammps-version.txt
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
    if [ "${base}" != "${base%.traj}" ] ; then
      # System force is now total force
      sed 's/fs_/ft_/g' < ${base} > ${TMPDIR}/${base}
      mv -f ${TMPDIR}/${base} ${base}
    fi
    ${SPIFF} -r 1e-${DIFF_PREC} $f $base > "$base.diff"
    RETVAL=$?
    if [ $RETVAL -ne 0 ]
    then
      if [ ${base##*\.} = 'out' ]
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
          ${SPIFF} -r 1e-${LOW_PREC} $f $base > /dev/null
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
