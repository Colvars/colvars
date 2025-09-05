#!/bin/bash
# -*- sh-basic-offset: 2; sh-indentation: 2; -*-

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

dir=$1

if [ -f ${dir}/disabled ] ; then
  echo "Directory ${dir} skipped."
  exit 0
fi

if [ -f ${dir}/skip_test.sh ]; then
  bash ${dir}/skip_test.sh
  if [ $? -eq 0 ]; then
    echo "Directory ${dir} skipped."
    exit 0
  fi
fi

echo -e "Entering $(${TPUT_BLUE})${dir}$(${TPUT_CLEAR}) ..."
cd $dir

if [ ! -d AutoDiff ] ; then
  echo "Directory ${dir} does not have AutoDiff/, skipped."
  exit 0
fi

# Precision requested to pass (negative powers of ten)
DIFF_PREC=6
DIFF_ABS_PREC=12
# Minimum precision to be tested
MIN_PREC=1

# Match only the trajectory file and force files
trajFileNamePattern="*.colvars.traj"
forceFileNamePattern="*_forces_*.dat"

SUCCESS=1
for f in "AutoDiff"/*
do
  base=`basename $f`
  if [[ ! ( $base = $trajFileNamePattern || $base = $forceFileNamePattern ) ]]
  then
    echo "Skip file $base"
    continue
  fi
  echo "Compare file $base"
  ${SPIFF} -r 1e-${DIFF_PREC} $f $base > "$base.diff"
  RETVAL=$?
  if [ $RETVAL -ne 0 ]
  then
    echo -e "\n*** Failure for file $(${TPUT_RED})$base$(${TPUT_CLEAR}): see `pwd`/$base.diff "
    SUCCESS=0
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
      SUCCESS=0
      echo " --> Passes at reduced precision 1e-${LOW_PREC}"
      break
    else
      # Test absolute error
      ${SPIFF} -m -a 1e-${DIFF_ABS_PREC} $f $base > /dev/null
      RETVAL=$?
      if [ $RETVAL -eq 0 ]
      then
        SUCCESS=1
        echo " --> Passes at absolute precision 1e-${DIFF_ABS_PREC}"
      else
        SUCCESS=0
        echo " --> Fails at minimum tested precision 1e-${LOW_PREC}"
        break
      fi
    fi
  else
    SUCCESS=1
  fi
done

if [ $SUCCESS -eq 1 ]
then
  echo "$(${TPUT_GREEN})Succeeded.$(${TPUT_CLEAR})"
  exit 0
else
  echo "$(${TPUT_RED})Failed.$(${TPUT_CLEAR})"
  exit 1
fi
