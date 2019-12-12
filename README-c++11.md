## Availability of Colvars features on different C++ versions

The majority of the Colvars module can be built with all major versions of the C++ language.  A few recent features are an exception, and require an executable built with C++11, which all recent compilers support but not all enable by default yet.

Because NAMD, LAMMPS or VMD do not support C++11 on all platforms yet, they are often built without those features.  This page points to the relevant information to identifying and potentially solving this issue.


### Specific Colvars features that require C++11

Currently the following variable types are only available when the code is built with C++11 standard or higher:
- `gspath` and `gzpath`
- `gspathCV` and `gzpathCV`
- `aspathCV` and `azpathCV`

Starting from 2019-06-02 `customFunction` also requires C++11, due to improvements in the Lepton library available from [the OpenMM repository](https://github.com/openmm/openmm).

### Status of C++ support in MD engines (as of 2019-10-17)

- _GROMACS_ already follows the C++11 standard in full.

- _LAMMPS_ is built by default with C++11, including most of the [precompiled builds](https://lammps.sandia.gov/download.html).  The C++11 transition is best handled through the [CMake build recipe](https://lammps.sandia.gov/doc/Build_cmake.html), which will automatically include or exclude C++11-specific sections of Colvars.  *Tip:* with CMake version 3.1 or later, try using the `-D CMAKE_CXX_STANDARD=11` flag to enable C++11 (if not already the compiler's default).  Build recipes following traditional `make` available in `src/MAKE/MACHINES` may require manual edits; check the README of the Colvars library in the LAMMPS distribution for details.

- _NAMD_ [precompiled builds](http://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=NAMD) will soon follow the C++11 standard but do not yet at the moment.  Build recipes for the Intel, Cray and IBM compilers in the `arch` folder are already set to require C++11.  Recipes based on `g++` should be adapted by replacing `-std=c++0x` with `-std=c++11`, to ensure that the `g++` version being used supports all needed features (in older GNU compilers, `0x` and `11` are not the same).

- _VMD_ currently does not provide build recipes with C++11 support.  Enabling C++11-dependent features will most likely require a custom build with manual edits to the `configure` script.
