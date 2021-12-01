#!/bin/bash
# -*- sh-basic-offset: 2; sh-indentation: 2; -*-

# Script to launch the regression tests
# It's best to have Gromacs compiled in double precision
# Reference files have been generated with Gromacs version 2020.3

gen_ref_output=''
TMPDIR=/tmp
DIRLIST=''
BINARY=gmx_d
# Default binary name when the option -DGMX_BUILD_MDRUN_ONLY=on is set during cmake.
# gmx_mpi_d should work too
BINARY_MPI=mdrun_mpi_d
SPIFF=spiff


while [ $# -ge 1 ]; do
  if { echo $1 | grep -q gmx ; }; then
    echo "Using GROMACS executable from $1"
    BINARY=$1
  elif [ "x$1" = 'x-g' ]; then
    gen_ref_output='yes'
    echo "Generating reference output"
  elif [ "x$1" = 'x-h' ]; then
    echo "Usage: ./run_tests.sh [-h] [-g] [path_to_gmx] [testdir1 [testdir2 ...]]"  >& 2
    echo "    The -g option (re)generates reference outputs in the given directories" >& 2
    echo "    If no executable is given, \"gmx_d\" is used" >& 2
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
else
  SPIFF=$(${TOPDIR}/devel-tools/get_spiff)
  if [ $? != 0 ] ; then
      echo "Error: could not be downloaded/built." >& 2
      echo "Using standard `spiff` command" >& 2
  else
      echo "Using spiff executable from $SPIFF"
      hash -p ${SPIFF} spiff
  fi
fi

if ! { echo ${DIRLIST} | grep -q 0 ; } then
  DIRLIST=`eval ls -d [0-9][0-9][0-9]_*`
fi

NUM_THREADS=3
NUM_CPUS=$(nproc)
if [ ${NUM_THREADS} -gt ${NUM_CPUS} ] ; then
  NUM_THREADS=${NUM_CPUS}
fi

TPUT_RED='true'
TPUT_GREEN='true'
TPUT_BLUE='true'
TPUT_CLEAR='true'
if which tput >& /dev/null ; then
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
  for script in test.tpr test.restart.tpr ; do
    for f in ${script%.tpr}.*diff; do if [ ! -s $f ]; then rm -f $f; fi; done # remove empty diffs only
    rm -f ${script%.tpr}.*{BAK,old,backup}
    for f in ${script%.tpr}.*{state,state.dat,state.stripped,out,traj,histogram?.dat,histogram?.dx,corrfunc.dat,pmf}
    do
      if [ ! -f "$f.diff" ]; then rm -f $f; fi # keep files that have a non-empty diff
    done
    rm -f *.xtc *.trr *.edr *.cpt *.gro *.log \#*${script%.tpr}.*  #Gromacs files
    rm -f
    rm -f *.out *.out.diff # Delete output files regardless
    rm -f *.ndx *.xyz
    rm -f test.dat
    if [ -f cleanup.sh ]
    then
      ./cleanup.sh
    fi
  done
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

  SCRIPTS="../Common/test.tpr ../Common/test.restart.tpr"
  restart=false

  cleanup_files



  # run simulation(s)
  for script in ${SCRIPTS} ; do

    basename=`basename ${script}`
    basename=${basename%.tpr}

    #prefix of the Gromacs generated files
    output=""

    #Different command line if it's a restart (checkpoint needed)
    if echo "${basename}" |grep -q "restart"
    then
      restart=true
    fi

    # Input files
    # Symbolink link to the colvars input file
    ln -sf ../Common/test.in test.dat

    # Handle MPI cases
    if [[ $dir  == *"tmpi"* ]]
    then
      COMMAND="${BINARY} mdrun"
      options="-ntmpi 2 -ntomp 2"
    elif [[ $dir == *"MPI"* ]]
    then
      COMMAND="mpirun -np 2 ${BINARY_MPI}"
      options="-ntomp 2"
    else # serial
      COMMAND="${BINARY} mdrun"
      options="-ntmpi 1 -ntomp 2"
    fi


    if [ -f "run.sh" ]
    then
      echo "Run special script."
      ./run.sh
      output=${basename}
    # Try running the test serial
    elif [ "$restart" = "false" ]
    then
      $COMMAND -s ${script} -deffnm ${basename} -colvars test.dat ${options} &> ${basename}.out
      RETVAL=$?
      output=${basename}
    else
      $COMMAND -s ${script} -deffnm ${basename} -noappend -cpi ${basename%.restart}.cpt -colvars test.dat -colvars_restart test.colvars.state.dat ${options} &> ${basename}.out
      RETVAL=$?
      output="${basename}.part0002"
    fi


    # Output of Colvars module, minus the version numbers
    grep "^colvars:" ${basename}.out | grep -v 'Initializing the collective variables module' \
      | grep -v 'Using GROMACS interface, version' > ${basename}.colvars.out


    if [ -f ${output}.colvars.state ] ; then
      # Filter out the version number from the state files to allow comparisons
      grep -sv 'version' ${output}.colvars.state > ${TMPDIR}/${output}.colvars.state.stripped
      mv -f ${TMPDIR}/${output}.colvars.state.stripped ${output}.colvars.state.stripped
      #Create symlink for restart input
      ln -sf ${output}.colvars.state{,.dat}
    fi

    #Convert the output of the traj files

    # If this test is used to generate the reference output files, copy them
    if [ "x${gen_ref_output}" = 'xyes' ]; then
      grep ':-) GROMACS -' ${basename}.out | head -n 1 > gromacs-version.txt
      grep 'Initializing the collective variables module, version' ${basename}.out | head -n 1 >> gromacs-version.txt
      grep 'Using GROMACS interface, version' ${basename}.out | head -n 1 >> gromacs-version.txt
      cp ${output}.colvars.state.stripped AutoDiff/
      cp ${output}.colvars.traj           AutoDiff/
      cp ${basename}.colvars.out            AutoDiff/
      if [ -f ${output}.histogram1.dat ] ; then
        cp -f ${output}.histogram1.dat AutoDiff/
      fi
      if [ -f ${output}.pmf ] ; then
        cp -f ${output}.pmf AutoDiff/
      fi
    fi

    if [ -f "run.sh" ]
    then
      break
    fi

  done

  # # now check results
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

  # # TODO: at this point, we may use the diff file to update the reference tests for harmless changes
  # # (e.g. keyword echos). Before then, figure out a way to strip the formatting characters produced by spiff.

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
