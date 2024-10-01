#!/bin/bash
# -*- sh-basic-offset: 2; sh-indentation: 2; -*-

# Script to launch the REMD regression tests
# It's best to have Gromacs compiled in double precision
# Reference files have been generated with Gromacs version 2024-dev

# if script called from `run_tests.sh` in the parent folder, retrieve the correct command binary from the argument
if [ "$1" ]
then
    BINARY=$1
else
    BINARY="gmx_mpi_d"
fi

mpirun -np 16 -oversubscribe $BINARY mdrun -multidir a b c d -s test.tpr -deffnm test -replex 10 -reseed 376 >& test.out

# Copy result files to be tested in the main folder
cp a/test.colvars.traj a/test.colvars.state .
cp a/test.log ./test.out

