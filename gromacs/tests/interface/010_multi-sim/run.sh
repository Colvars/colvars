#!/bin/bash
# -*- sh-basic-offset: 2; sh-indentation: 2; -*-

# Script to launch the REMD regression tests
# It's best to have Gromacs compiled in double precision
# Reference files have been generated with Gromacs version 2020.6

BINARY_MPI=mdrun_mpi_d

# Attempt replica exchange every ten steps
mpirun  -np 16 $BINARY_MPI -multidir a b c d -s test.tpr  -deffnm test -replex 10 -reseed 376 -colvars test.dat &> test.out
# Copy result files to be tested in the main folder
cp a/test.colvars.traj a/test.colvars.state .
