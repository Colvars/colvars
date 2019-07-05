## Availability of Colvars features on different C++ versions

The majority of the Colvars module can be built with all major versions of the C++ language.  A few recent features are an exception, and require an executable built with C++11, which all recent compilers support but none enable by default yet.

Because NAMD, LAMMPS or VMD do not support C++11 on all platforms yet, they are often built without those features.  This page points to the relevant information to identifying and potentially solving this issue.


### Specific Colvars features that require C++11

Currently, the `gspath`, `gzpath`, `gspathCV` and `gzpathCV` collective variables are only available when the code is built with C++11 standard or higher.


### Status of C++ support in MD engines (as of 2019-07-02)

- _GROMACS_ already follows the C++11 standard in full.

- _LAMMPS_ can be built for the most part with C++11, and some of the [precompiled builds](https://lammps.sandia.gov/download.html) already support it.  Current build recipes following traditional `make` available in `src/MAKE/MACHINES` require manual edits.  It is often easier to use the [CMake build recipe](https://lammps.sandia.gov/doc/Build_cmake.html), and specify manually the `-std=c++11` compiler flag or simply enable one of the optional packages that require C++11 (e.g. `KOKKOS`).

- _NAMD_ [precompiled builds](http://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=NAMD) will soon follow the C++11 standard but do not yet at the moment.  Build recipes for the Intel, Cray and IBM compilers in the `arch` folder are already set to require C++11.  Recipes based on `g++` should be adapted by replacing `-std=c++0x` with `-std=c++11` to ensure that the `g++` version being used supports all needed features.

- _VMD_ currently does not provide build recipes with C++11 support.  Enabling C++11-dependent features will most likely require a custom build with manual edits to the `configure` script.
