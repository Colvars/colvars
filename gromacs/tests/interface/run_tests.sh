#!/bin/bash

TOPDIR=$(git rev-parse --show-toplevel)
if [ ! -d ${TOPDIR} ] ; then
    echo "Error: cannot identify top project directory." >& 2
    exit 1
fi

MPI_TESTS=[0-9][0-9][0-9]_MPI_*

tests=([0-9][0-9][0-9]_*)

for DIR in $MPI_TESTS
do
  tests=(${tests[@]/$DIR})
done

ALL_SUCCESS=1

echo "Running tests ${tests[@]}"

../library/run_tests.sh $1 ${tests[@]}

if [ $? -ne 0 ]
then
    ALL_SUCCESS=0
fi

# Run tests that depend on an MPI build
if source ${TOPDIR}/devel-tools/load-openmpi.sh ; then
    for DIR in $MPI_TESTS
    do
      if pushd $DIR > /dev/null ; then
          ./run.sh $1
          if [ $? -ne 0 ]
          then
              ALL_SUCCESS=0
          fi
          popd > /dev/null
      fi
    done
fi


if [ $ALL_SUCCESS -eq 1 ]
then
  echo "$(${TPUT_GREEN})All tests succeeded.$(${TPUT_CLEAR})"
  exit 0
else
  echo "$(${TPUT_RED})There were failed tests.$(${TPUT_CLEAR})"
  exit 1
fi
