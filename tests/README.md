## Colvars library tests

This folder contains example configuration files that define computations via the Colvars library.  Most configurations are written using features that are supported by multiple backends (GROMACS, LAMMPS, NAMD and VMD).

Because not all engines support name-based selections or PDB file-based selections, and in the interest of compactness and portability, the primary atom selection method used here is index groups, which will require one or more `indexFile` keywords to be included.  For other selection methods, [engine-specific tests](#Engine-specific-tests) may be also created without using the tools from this folder.


### Test input syntax and dependencies

Only a Colvars configuration file must be created: in each engine's test folder, there are already a topology for a small molecular system (currently, deca-alanine in a vacuum) and an index file named `index.ndx`.

The configuration files included in this folder use index groups with one or more of these:
- Groups named `group1`, `group2` ... `group10`, which are suitable when a variable is defined from multiple groups; each contains 4 atoms or more.
- The `RMSD_atoms` group, which is suitable for variables that require a set of reference coordinates (either directly or indirectly, e.g. through `fittingGroup` and similar).
- The `heavy_atoms` group, which is a suitable choice for most cases.  Note that at least in biomolecular systems, leaving out hydrogen atoms from a variable's definition minimizes the chance of tests failing due to large fluctuations.

For the above two groups, please use the file names `rmsd_atoms_refpos.xyz` and `rmsd_atoms_refpos.xyz` to supply reference coordinates, since these files are already set up in each engine's test folder.


### Testing procedure and pass/fail criteria

For each backend, there is a folder `<engine>/tests/library`, which contains a `run_tests.sh` script and full input decks for each of the tests included here.  Numerical results are compared by that script against previously recorded values (i.e. regression tests are performed).

Each test runs a short MD simulation (e.g. 20 steps) with the program exiting gracefully, followed by a second MD run of the same length, restarted from the first.  Colvars trajectory files, state files, and other output files are all compared against reference values: if a relative deviation greater than 1.0e-6 is observed, that test fails.

Tests are run in double precision wherever this is feasible: therefore, most tests match their expected results within a far smaller error than 1.0e-6.  However, remember that a *full MD simulation* is performed for each test, and changes in the MD engine's version and choices of compiler may introduce differences whichever version of Colvars is used.

In the automated tests performed on GitHub, a CentOS 7 Linux environment and the GCC 4.8 compiler are used to build Colvars and other codes (exception: GROMACS uses GCC 10 from the same environment).  A text file in each test folder states which code versions were used to generate the reference output.  Occasionally, reference outputs may need to be updated to reflect changes in codes unrelated to Colvars: this is done as rarely and conservatively as possible.

Please keep the above in mind if you plan on running tests in your own environment for development purposes, and feel free to reach out via GitHub issue if a test's failure is suspicious.


### Engine-specific tests

Some Colvars features are specific to a certain engine: for example, not all codes implement atom name-based or PDB file-based selections, or Colvars may be using code from that engine to perform a computation.  For those cases, additional tests are also found in the `<engine>/tests/library` folder.  There is also a `<engine>/tests/interface` folder to test functionality that is unique to that code.

*Aside from the above exceptions, it is highly recommended to create test inputs at the top level, i.e. inside this folder, so that they may be reused and their features tested in every scenario.*


### Adding a new test

- First, a Colvars configuration file should be created; the default extension used is `.in` for simplicity (but any extension would otherwise work for a Colvars file).  The file may either be complete on its own:
```
nano my_config_name.in
```
or be split up between a colvar and a bias segment:
```
nano colvar_name.in
nano bias_name.in
```
In either case, a file called simply `test.in` will be generated.  This file will be placed inside a folder dedicated to running the corresponding test.  This folder will contain any additional input files, as well as the reference output files for the test.

- Next, add the following lines to the script `build_tests.sh` to create a test folder:
```
create_test_dir <my_config_name>
```
where for a colvar-bias file scheme it would be a good idea to simply use `<colvar_name>_<bias_name>`.  Immediately after, add this line to create the input file in it, using either the single-file scheme:
```
write_colvars_config <my_config_name> ""
```
or the colvar-bias scheme:
```
write_colvars_config <colvar_name> <bias_name>
```
(If the new test has additional dependencies other than the Colvars config file, those would have to be copied manually as well.)

- Lastly, run the script to execute these tasks:
```
./build_tests.sh ../<engine>/tests/library
```
and navigate to the test's parent folder to run it for the first time and generate reference files:
```
cd ../<engine>/tests/library
./run_tests.sh -g <path_to_executable> 000_<my_config_name>
```

- It is always a good idea to check that the test matches the reference files just generated:
```
./run_tests.sh <path_to_executable> 000_<my_config_name>
```
