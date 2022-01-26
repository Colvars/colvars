## Availability of Colvars features on different C++ versions

The majority of the Colvars module can be built with all major versions of the C++ language.  A few recent features are an exception, and require an executable built with C++11, which all recent compilers support but not all enable by default yet.

At this time, the issue only affects VMD, whose build system aims to keep supporting also some Linux versions that are no longer under active support.

Upon initialization, Colvars will always print the version of C++ with which it was built, so that you may compare it with the information in this README.


### Specific Colvars features that require C++11

Currently the following variable types are only available when the code is built with C++11 standard or higher:
- `gspath` and `gzpath`
- `gspathCV` and `gzpathCV`
- `aspathCV` and `azpathCV`

Starting from 2019-06-02 `customFunction` also requires C++11, due to improvements in the Lepton library available from [the OpenMM repository](https://github.com/openmm/openmm).

### Status of C++ support in MD engines (as of 2021-01-21)

- **GROMACS** currently follows the C++14 standard, which is backward compatible with C++11.

- **LAMMPS** follows the C++11 standard.

- **NAMD** follows the C++11 standard, required by Charm++.

- **VMD** precompiled builds *do not* use C++11 yet.  If you are building your own VMD, try using the optional `CXX11` flag when running the `configure` script, which will activate C++11 flags for the `LINUXAMD64` architecture: however, support for this option is still a work in progress.  Alternatively, you may use the custom `c++11` patch from [this repository](https://github.com/giacomofiorin/vmd-patches/) as a temporary fix.
