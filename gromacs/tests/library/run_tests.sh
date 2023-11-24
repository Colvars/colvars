#!/bin/bash
# -*- sh-basic-offset: 2; sh-indentation: 2; -*-

# Script to launch the regression tests
# It's best to have Gromacs compiled in double precsision
# Reference files have been generated with Gromacs version 2020.3

gen_ref_output=''
TMPDIR=/tmp
DIRLIST=''
BINARY=gmx_d
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
  for script in test.tpr test.restart.tpr ; do
    for f in ${script%.tpr}.*diff; do if [ ! -s $f ]; then rm -f $f; fi; done # remove empty diffs only
    rm -f ${script%.tpr}.*{BAK,old,backup}
    for f in ${script%.tpr}.*{state,state.dat,state.stripped,out,traj,histogram?.dat,histogram?.dx,corrfunc.dat,pmf}
    do
      if [ ! -f "$f.diff" ]; then rm -f $f; fi # keep files that have a non-empty diff
    done
    rm -f *.xtc *.trr *.edr *.cpt *.gro *.log \#*${script%.tpr}.*  #Gromacs files
    rm -f
    rm -f *.out *.out.diff *.err # Delete output files regardless
    rm -f *.ndx *.xyz
    rm -f test.dat
  done
}


if ${BINARY} mdrun -colvars 2>&1 | grep -q 'Unknown command-line option -colvars' ; then
  MDMODULES_INTERFACE=yes
else
  MDMODULES_INTERFACE=no
fi

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

  # Run simulation(s)
  for basename in test test.restart ; do

    # Input files
    # Symbolink link to the colvars input file, index file, and xyz file
    ln -sf test.in test.dat
    ln -sf ../Common/da.ndx index.ndx
    if [ -f test.dat ] && grep -q "refPositionsFile rmsd_" test.dat
    then
        ln -fs ../Common/rmsd_atoms_refpos.xyz ./
        ln -fs ../Common/rmsd_atoms_refpos2.xyz ./
        ln -fs ../Common/rmsd_atoms_random.xyz ./
    fi
    if [ -f test.dat ] && grep -q "refPositionsFile heavy_" test.dat
    then
        ln -fs ../Common/heavy_atoms_refpos.xyz heavy_atoms_refpos.xyz
    fi

    # Try running the test

    if [ "x${MDMODULES_INTERFACE}" == "xyes" ] ; then

      if [ "${basename}" == "test" ] ; then
        ${BINARY} grompp -f ../Common/test.mdp -c ../Common/da.pdb -p ../Common/da.top -t ../Common/da.trr -o ${basename}.tpr 2> ${basename}.grompp.err 1> ${basename}.grompp.out
        ${MPIRUN_CMD} ${BINARY} mdrun ${TMPI_TASKS} -s ${basename}.tpr -ntomp ${NUM_THREADS} -deffnm ${basename} 2> ${basename}.err 1> ${basename}.out
        RETVAL=$?
      fi

      if [ "${basename}" == "test.restart" ] ; then

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
          ${BINARY} grompp -f ${NEW_MDP} -c ../Common/da.pdb -p ../Common/da.top -t ${basename%.restart}.cpt -o ${basename}.tpr 2> ${basename}.grompp.err 1> ${basename}.grompp.out
          rm -f ${NEW_MDP} ${NEW_CVCONF}
          ${MPIRUN_CMD} ${BINARY} mdrun ${TMPI_TASKS} -s ${basename}.tpr -ntomp ${NUM_THREADS} -deffnm ${basename} -noappend 2> ${basename}.err 1> ${basename}.out
          RETVAL=$?

          output=${basename}.part0001
          for file in ${output}.* ; do
            # Remove the part number
            mv -f ${file} ${file/.part0001/}
          done

        else

          # Restart both GROMACS and Colvars using the GROMACS checkpoint file
          ${BINARY} convert-tpr -s ${basename%.restart}.tpr -nsteps 40 -o ${basename}.tpr 2> ${basename}.grompp.err 1> ${basename}.grompp.out
          ${MPIRUN_CMD} ${BINARY} mdrun ${TMPI_TASKS} -s ${basename}.tpr -ntomp ${NUM_THREADS} -deffnm ${basename} -noappend -cpi ${basename%.restart}.cpt 2> ${basename}.err 1> ${basename}.out

          RETVAL=$?
          output=${basename}.part0002
          for file in ${output}.* ; do
            # Remove the part number
            mv -f ${file} ${file/.part0002/}
          done

        fi
      fi

    else

      if [ "${basename}" == "test" ] ; then
        ln -fs ${basename}.in ${basename}.dat
        ${MPIRUN_CMD} ${BINARY} mdrun ${TMPI_TASKS} -s ../Common/${basename} -deffnm ${basename} -colvars test.dat &> ${basename}.out
        RETVAL=$?
        ln -fs ${basename}.colvars.state{,.dat}
      fi

      if [ "${basename}" == "test.restart" ] ; then
        ${MPIRUN_CMD} ${BINARY} mdrun ${TMPI_TASKS} -s ../Common/${basename} -deffnm ${basename} -noappend -cpi ${basename%.restart}.cpt -colvars ${basename%.restart}.dat -colvars_restart ${basename%.restart}.colvars.state.dat &> ${basename}.out
        RETVAL=$?
        output="${basename}.part0002"
        for file in ${output}.* ; do
          # Remove the part number
          mv -f ${file} ${file/.part0002/}
        done
      fi

    fi

    # Filter out the version numbers to allow comparisons
    grep "^colvars:" ${basename}.log \
      | grep -v 'Initializing the collective variables module' \
      | grep -v 'Using GROMACS interface, version' > ${basename}.colvars.out
    if [ -f ${basename}.colvars.state ] ; then
      grep -sv 'version' ${basename}.colvars.state \
           > ${TMPDIR}/${basename}.colvars.state.stripped && \
        mv -f ${TMPDIR}/${basename}.colvars.state.stripped ${basename}.colvars.state.stripped
    fi

    # If this test is used to generate the reference output files, copy them
    if [ "x${gen_ref_output}" = 'xyes' ]; then
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

  # # now check results
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
