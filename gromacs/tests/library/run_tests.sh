#!/bin/bash
# -*- sh-basic-offset: 2; sh-indentation: 2; -*-

# Script to launch the regression tests
# It's best to have Gromacs compiled in double precsision
# Reference files have been generated with Gromacs version 2020.3

gen_ref_output=''
TMPDIR=/tmp
DIRLIST=''
BINARY=gmx_d


while [ $# -ge 1 ]; do
  if { echo $1 | grep -q gmx ; }; then
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
    # Symbolink link to the colvars input file, index file, and xyz file
    ln -sf test.in test.dat
    ln -sf ../Common/da.ndx index.ndx
    if grep -q "refPositionsFile rmsd_" test.dat
    then
        ln -fs ../Common/rmsd_atoms_refpos.xyz rmsd_atoms_refpos.xyz
    fi
    if grep -q "refPositionsFile heavy_" test.dat
    then
        ln -fs ../Common/heavy_atoms_refpos.xyz heavy_atoms_refpos.xyz
    fi

    # Try running the test
    if [ "$restart" = "false" ]
    then
      $BINARY mdrun -s ${script} -deffnm ${basename} -colvars test.dat &> ${basename}.out
      RETVAL=$?
      output=${basename}
    else
      if [ -f test.restart.in ] ; then
        ln -fs test.restart.in test.dat
      fi
      $BINARY mdrun -s ${script} -deffnm ${basename} -noappend -cpi ${basename%.restart}.cpt -colvars test.dat -colvars_restart test.colvars.state.dat &> ${basename}.out
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
      spiff -r 1e-${DIFF_PREC} $f $base > "$base.diff"
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
