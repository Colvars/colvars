## Gromacs mwABF regression test

This directory contains a test for multiple-walker ABF in Colvars/Gromacs.


### Version tested

The tests are valid against Gromacs version 2024.3 and following.


### Files tested

Only the traj files of the simulation **a** is being tested.


### Scripts

  - `create_tpr.sh` is used to generate the .tpr
  - `run.sh` is called by a `run_tests.sh` script to run the REMD simulation
  - `cleanup.sh` is called by a `run_tests.sh` to clean the simulations files in the subfolders.
