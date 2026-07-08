#!/bin/bash
# -*- sh-basic-offset: 2; sh-indentation: 2; -*-

# Script to launch the regression tests
# It's best to have Gromacs compiled in double precsision
# Reference files have been generated with Gromacs version 2020.3

TMPDIR=/tmp
DIRLIST=''
BINARY=gmx_d
SPIFF=spiff

cleanup='true'
verbose='false'
gen_ref_output='false'

while [ $# -ge 1 ]; do
  if { echo $1 | grep -q gmx ; }; then
    BINARY=$1
  elif [ "x$1" = 'x-g' ]; then
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
    exit 0
  else
    DIRLIST=`echo ${DIRLIST} $1`
  fi
  shift
done

echo "Using GROMACS executable: ${BINARY}"
if $verbose ; then
  echo "Will run tests in the following directories: ${DIRLIST}"
fi

TOPDIR=$(git rev-parse --show-toplevel)
if [ ! -d ${TOPDIR} ] ; then
  echo "Error: cannot identify top project directory." >& 2
else
  SPIFF=$(${TOPDIR}/devel-tools/get_spiff)
  if [ $? != 0 ] ; then
      echo "Error: could not be downloaded/built." >& 2
      echo "Using standard `spiff` command" >& 2
  else
      echo "Using spiff executable: $SPIFF"
      hash -p ${SPIFF} spiff
  fi
fi

if ! { echo ${DIRLIST} | grep -q 0 ; } then
  DIRLIST=`eval ls -d [0-9][0-9][0-9]_*`
fi


NUM_TASKS=${NUM_TASKS:-4}
NUM_THREADS=${NUM_THREADS:-2}

MPI_BUILD=no
TMPI_TASKS="-ntmpi ${NUM_TASKS}"
MPIRUN_CMD=""
if echo $BINARY | grep -q gmx_mpi ; then
  MPI_BUILD=yes
  if source ${TOPDIR}/devel-tools/load-openmpi.sh ; then
    MPIRUN_CMD="mpirun -n ${NUM_TASKS} -oversubscribe"
    echo "Will run mdrun as: ${MPIRUN_CMD} ${BINARY} mdrun -ntomp ${NUM_THREADS}"
    TMPI_TASKS=""
  fi
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
DIFF_PREC=6
# Minimum precision to be tested
MIN_PREC=1

cleanup_files() {
  local -a scripts=("$@")
  if (( $# == 0 )); then
    scripts=(*.tpr)
  fi
  for script in "${scripts[@]}"; do
    base="${script%.tpr}"
    for f in ${base}.*diff; do if [ ! -s $f ]; then rm -f $f; fi; done # remove empty diffs only
    rm -f ${base}.*{BAK,old,backup}
    for f in ${base}.*{state,state.dat,state.stripped,out,traj,histogram?.dat,histogram?.dx,corrfunc.dat,pmf}
    do
      if [ ! -f "$f.diff" ]; then rm -f $f; fi # keep files that have a non-empty diff
    done
    rm -f *.xtc *.trr *.edr *.cpt *.gro *.log mdout.mdp \#*  #Gromacs files
    rm -f metadynamics1.*.txt replicas.registry.txt *.hills *.tpr
    rm -f *.out *.out.diff *.err # Delete output files regardless
    rm -f *.ndx *.xyz
    rm -f ${base}.dat
  done
  # if cleanup.sh exists, run it to clean up any additional files
  if [ -f cleanup.sh ]; then
    ./cleanup.sh
  fi
}


for dir in ${DIRLIST} ; do

  dir="${dir%/}" # Remove trailing / if present
  if [ -f ${dir}/disabled ] ; then
    if $verbose ; then
      echo "Skipping disabled test $dir"
    fi
    continue
  fi

  # Skip MPI tests if the binary is not an MPI build
  if [ $MPI_BUILD = no ] && echo ${dir} | grep -q MPI ; then
    echo "Skipping $(${TPUT_BLUE})${dir}$(${TPUT_CLEAR}) (MPI test, but not an MPI build)"
    continue
  fi

  echo -ne "Entering $(${TPUT_BLUE})${dir}$(${TPUT_CLEAR}) ... "
  cd $dir

  if [ ! -d AutoDiff ] ; then
    echo ""
    echo "  Creating directory AutoDiff, use -g to fill it."
    mkdir AutoDiff
    cd $BASEDIR
    continue
  else

    if $gen_ref_output; then

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

  if $cleanup; then
    cleanup_files
  fi

  simulations=(test test.restart)
  if [ "${dir##*/}" = "000_multiple_walkers_mtd" ] ; then
    simulations=(test{.rep1,.rep2} test{.rep1,.rep2}.restart)
  fi

  # Run simulation(s)
  if [ -f run.sh ]
  then
    # Special run script e.g. for interface tests
    ./run.sh $BINARY
    RETVAL=$?

    # For multi-replica simulations, copy output from first replica
    if [ -d a ] && echo ${dir} | grep -q MPI; then
      basename=test
      logfile=a/${basename}.log
      # Filter out the version numbers to allow comparisons
      grep "^colvars:" ${logfile} \
        | grep -v 'Initializing the collective variables module' \
        | grep -v 'Using GROMACS interface, version' > ${basename}.colvars.out

      if [ -s a/${basename}.colvars.state ] ; then
        grep -sv 'version' a/${basename}.colvars.state \
            > ${TMPDIR}/${basename}.colvars.state.stripped && \
          mv -f ${TMPDIR}/${basename}.colvars.state.stripped ${basename}.colvars.state.stripped
      fi

      cp a/${basename}.colvars.traj .
    fi
  else
    for basename in ${simulations[@]} ; do

      # Create symlinks to the Colvars config file, index file, and xyz file
      if [ -f ${basename}.in ] ; then
        ln -sf ${basename}.in ${basename}.dat
        ln -sf ../Common/da.ndx index.ndx
        if grep -q "refPositionsFile rmsd_" ${basename}.dat
        then
          ln -fs ../Common/rmsd_atoms_refpos.xyz ./
          ln -fs ../Common/rmsd_atoms_refpos2.xyz ./
          ln -fs ../Common/rmsd_atoms_random.xyz ./
        fi
        if grep -q "refPositionsFile heavy_" ${basename}.dat
        then
          ln -fs ../Common/heavy_atoms_refpos.xyz heavy_atoms_refpos.xyz
        fi
      fi

      # Try running the test

      structure_file="../Common/da.pdb"
      topology_file="../Common/da.top"
      initial_trajectory_file="../Common/da.trr"
      if [ -f get_structure_file.sh ] ; then
        structure_file=$(./get_structure_file.sh)
      fi
      if [ -f get_topology_file.sh ] ; then
        topology_file=$(./get_topology_file.sh)
      fi
      if [ -f get_initial_trajectory_file.sh ] ; then
        initial_trajectory_file=$(./get_initial_trajectory_file.sh)
      fi
      if [ "${basename%.restart}" == "${basename}" ] ; then
        # Initial run
        MDP=../Common/test.mdp
        if [ -f ${basename}.mdp ] ; then
          MDP=${basename}.mdp
        fi
        if [[ -z "$initial_trajectory_file" ]] ; then
          ${BINARY} grompp -f ${MDP} -c $structure_file -p $topology_file -o ${basename}.tpr 2> ${basename}.grompp.err 1> ${basename}.grompp.out
        else
          ${BINARY} grompp -f ${MDP} -c $structure_file -p $topology_file -t $initial_trajectory_file -o ${basename}.tpr 2> ${basename}.grompp.err 1> ${basename}.grompp.out
        fi
        ${MPIRUN_CMD} ${BINARY} mdrun ${TMPI_TASKS} -s ${basename}.tpr -ntomp ${NUM_THREADS} -deffnm ${basename} 2> ${basename}.err 1> ${basename}.out
        RETVAL=$?
      fi

      if [ "${basename%.restart}" != "${basename}" ] ; then
        # Restart run

        if [ -n "${FORCE_INPUT_STATE_FILE}" ] ; then

          # Restart GROMACS using the checkpoint but Colvars using its own state file

          # Add defaultInputStateFile to the Colvars config
          NEW_CVCONF=$(mktemp test.XXXXX.in)
          cat test.in > ${NEW_CVCONF}
          echo "defaultInputStateFile test.colvars.state" >> ${NEW_CVCONF}

          NEW_MDP=$(mktemp test.XXXXX.mdp)
          cat ../Common/test.mdp > ${NEW_MDP}
          sed -i "s/test.in/${NEW_CVCONF}/" ${NEW_MDP}
          # Mimic the initial step of a job restarted from checkpoint, to be
          # consistent with reference outputs
          echo "init-step = 20" >> ${NEW_MDP}
          ${BINARY} grompp -f ${NEW_MDP} -c $structure_file -p $topology_file -t ${basename%.restart}.cpt -o ${basename}.tpr 2> ${basename}.grompp.err 1> ${basename}.grompp.out
          rm -f ${NEW_MDP} ${NEW_CVCONF}
          ${MPIRUN_CMD} ${BINARY} mdrun ${TMPI_TASKS} -s ${basename}.tpr -ntomp ${NUM_THREADS} -deffnm ${basename} -noappend 2> ${basename}.err 1> ${basename}.out
          RETVAL=$?
        else

          # Restart both GROMACS and Colvars using the GROMACS checkpoint file
          ${BINARY} convert-tpr -s ${basename%.restart}.tpr -nsteps 40 -o ${basename}.tpr 2> ${basename}.grompp.err 1> ${basename}.grompp.out
          ${MPIRUN_CMD} ${BINARY} mdrun ${TMPI_TASKS} -s ${basename}.tpr -ntomp ${NUM_THREADS} -deffnm ${basename} -noappend -cpi ${basename%.restart}.cpt 2> ${basename}.err 1> ${basename}.out
          RETVAL=$?
        fi
      fi

      logfile=${basename}.log
      if [ -s ${basename}.part0002.log ] ; then
        logfile=${basename}.part0002.log
      fi

      # Filter out the version numbers to allow comparisons
      grep "^colvars:" ${logfile} \
        | grep -v 'Initializing the collective variables module' \
        | grep -v 'Using GROMACS interface, version' > ${basename}.colvars.out

      if [ -s ${logfile%.log}.colvars.state ] ; then
        grep -sv 'version' ${logfile%.log}.colvars.state \
            > ${TMPDIR}/${logfile%.log}.colvars.state.stripped && \
          mv -f ${TMPDIR}/${logfile%.log}.colvars.state.stripped ${logfile%.log}.colvars.state.stripped
      fi

      # If this test is used to generate the reference output files, copy them
      if $gen_ref_output; then
        grep ':-) GROMACS -' ${basename}.out | head -n 1 > gromacs-version.txt
        grep 'Initializing the collective variables module, version' ${basename}.out | head -n 1 >> gromacs-version.txt
        grep 'Using GROMACS interface, version' ${basename}.out | head -n 1 >> gromacs-version.txt
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

    done
  fi

  for file in * ; do
    # Remove part numbers
    if [ ${file} != ${file/.part0001/} ] ; then
      mv -f ${file} ${file/.part0001/}
    fi
    if [ ${file} != ${file/.part0002/} ] ; then
      mv -f ${file} ${file/.part0002/}
    fi
  done

  # # now check results
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
      cleanup_files ${simulations[@]}
    fi
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
