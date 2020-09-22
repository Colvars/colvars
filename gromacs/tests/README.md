## Gromacs regression tests

This folder contains regression tests input and output files for Gromacs, testing various user-facing features of Colvars.
The tests under `library` have been generetad to match at best the NAMD ones. The system and simulation parameters have been chosen to reflect NAMD ones (SHAKE algorithm, time step of 1ps, use of chosen velocities). PBC use is mandatory (unlike NAMD) due to the Verlet algorithm. See `library/Common/create_tpr.sh` for more information.

Specific tests of the Colvars-Gromacs interface are located in `interface`. Here, simulations parameters are different to be able to test MPI parallelization.

### Version tested

The tests are valid against Gromacs version 2020.X. The inputs (.tpr) and reference files have been generated with 2020.3.

### Default usage

The default task of the `run_tests.sh` script (found in `library` and `interface`) is to run the input files contained in each directory, and compare their result to the files in the corresponding `AutoDiff` subfolder.  The latter must have been created previously using `./run_tests.sh -g`.

### Script usage

This is the usage information for `run_tests.sh`:
```
Usage: ./run_tests.sh [-h] [-g] [path_to_namd2] [testdir1 [testdir2 ...]]
    The -g option (re)generates reference outputs in the given directories
    If no executable is given, "namd2" is used
    If no directories are given, all matches of [0-9][0-9][0-9]_* are used
```

### File organization

All tests matching `000_*` use the same Gromacs input scripts `Common/test.mdp` and `Common/test.restart.mdp`.

Any input files that are meant to be reused should be also located in `Common`.

### Known issues

Floating-point trajectories are obviously sensitive to the compiler version and optimization flags used, and their deviations may exceed the tolerance normally used (1.0e-6).  For this reason, `run_tests.sh` also prints results of the tests at higher tolerance level.
