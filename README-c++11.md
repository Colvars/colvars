## Availability of Colvars features on different C++ versions

The majority of the Colvars module can be built with all major versions of the C++ language.  A few recent features are an exception, and require an executable built with C++11, which all recent compilers support but not all enable by default yet.

At this time, the issue only affects VMD.  Because the VMD build system does not support C++11 on all platforms yet, VMD is often built without those features.  This page points to the relevant information to identifying and potentially solving this issue.


### Specific Colvars features that require C++11

Currently the following variable types are only available when the code is built with C++11 standard or higher:
- `gspath` and `gzpath`
- `gspathCV` and `gzpathCV`
- `aspathCV` and `azpathCV`

Starting from 2019-06-02 `customFunction` also requires C++11, due to improvements in the Lepton library available from [the OpenMM repository](https://github.com/openmm/openmm).

### Status of C++ support in MD engines (as of 2020-10-29)

- *GROMACS* currently follows the C++14 standard, which is backward compatible with C++11.

- *LAMMPS* follows the C++11 standard.

- *NAMD* [precompiled builds](http://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=NAMD) and most of the build recipes contained in the `arch` folder follow the C++11 standard.

- *VMD* currently does not provide precompiled builds or build recipes with C++11 support.  Enabling C++11-dependent features requires a custom build with a modified the `configure` script.  [This repository](https://github.com/giacomofiorin/vmd-patches) contains, among others, an example patch for making such changes.
