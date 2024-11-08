#!/bin/bash
# -*- sh-basic-offset: 2; sh-indentation: 2; -*-

# Script to launch the REMD regression tests
# It's best to have Gromacs compiled in double precision
# Reference files have been generated with Gromacs version 2024-dev

echo -n "Running $(basename $PWD)..."

#set -e

# if script called from `run_tests.sh` in the parent folder, retrieve the correct command binary from the argument
if [ -x "$1" ]
then
    BINARY=$1
else
    BINARY="gmx_mpi_d"
fi

for i in a b c d
do
  cd $i
  $BINARY grompp -f system.mdp -p ../../Common/system.top -c ../../Common/system.gro -o test.tpr &> grompp.log
  cd ..
done

mpirun -np 16 -oversubscribe $BINARY mdrun -multidir a b c d -s test.tpr -deffnm test -reseed 376 >& test.out

labels=(a b c d)
for number in $(seq 1 ${#labels[@]}) ; do
  index=$((${number} - 1))
  label=${labels[${index}]}
  replica_msg=$(grep -s "colvars: Enabling multiple replicas: this is replica number" ${label}/test.log)
  if [ -z "${replica_msg}" ] ; then
    echo "Error: missing replica initialization message." >& 2
    exit 1
  fi
  replica_num=$(echo ${replica_msg} | cut -d' ' -f 9)
  replicas_count=$(echo ${replica_msg} | cut -d' ' -f 11)
  replicas_count=${replicas_count%.}
  if [ ${replicas_count} != 4 ] || [ ${replica_num} != ${number} ] ; then
    echo -e "Error: wrong initialization for replica number ${number}; got the following message:\n\n${replica_msg}\n" >& 2
    exit 1
  fi
done

# Copy result files to be tested in the main folder
cp a/test.colvars.traj a/test.colvars.state .
cp a/test.log ./test.out
cp a/test.all.czar.pmf ./test.all.czar.pmf

echo " Done."

