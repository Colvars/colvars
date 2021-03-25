## NAMD regression tests

This folder contains regression tests input and output files for NAMD, testing various user-facing features of Colvars.  Specific tests of the Colvars-NAMD interface are located in `namd/tests/interface`.

### Default usage

The default task of the `run_tests.sh` script is to run the input files contained in each directory, and compare their result to the files in the corresponding `AutoDiff` subfolder.  The latter must have been created previously using `./run_tests.sh -g`.

### Script usage

This is the usage information for `run_tests.sh`:
```
Usage: ./run_tests.sh [-h] [-g] [path_to_namd2] [testdir1 [testdir2 ...]]
    The -g option (re)generates reference outputs in the given directories
    If no executable is given, "namd2" is used
    If no directories are given, all matches of [0-9][0-9][0-9]_* are used
```

### File organization

With very few exceptions, the Colvars configuration file in each folder is called `test.in`.

Occasionally, a file `test.legacy.in` is also available: this includes syntax that is now deprecated, but is kept around so that the same test could be repeated with much earlier code versions.

All tests prefixed by `000_` re-use the same NAMD input scripts `Common/test.namd` and `Common/test.restart.namd` and different versions of `test.in` containing features that are NAMD-agnostic (i.e. they are also used in LAMMPS and GROMACS).  The original Colvars configuration files are located at `../../../tests`.

Tests numbered sequentially starting from `001` use specific NAMD scripts.

Any input files that are meant to be reused in multiple tests should be located in `Common`.

### Known issues

- Floating-point trajectories are obviously sensitive to the compiler version and optimization flags used, and their deviations may exceed the tolerance normally used (1.0e-6).  For this reason, `run_tests.sh` also prints results of the tests at higher tolerance level.  Currently GCC 4.8 is used to build code to perform numerical tests at the highest precision.
