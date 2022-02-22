## Gromacs REMD regression test

This folder contains regression tests input and output files for a REMD (by temperature) simulation.
This REMD has specific conditions to test the correct exchange and that PBC are correctly treated:
  - 3 of the simulations (*b*, *c*, *d*) has the same starting point
  - the 4th (*a*) has a different starting point with one of the colvars atoms has jumped in the over side of the box.
  - for 2 of them (*a*, *b*), an harmonic constraint has been set between the 2 colvars atoms
  - for the others 2 (*c*, *d*), there isn't an harmonic constraint.


### Version tested

The tests are valid against Gromacs version 2020.X. The inputs (.tpr) and reference files have been generated with 2020.6.


### Files tested

Only the traj files of the simulation **a** is being tested.


### Scripts

  - `create_tpr.sh` is used to generate the .tpr
  - `run.sh` is called by a `run_tests.sh` script to run the REMD simulation
  - `cleanup.sh` is called by a `run_tests.sh` to clean the simulations files in the subfolders.
