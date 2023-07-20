## Availability of Colvars features depending on C++ language versions

The majority of the Colvars module can be built with all major versions of the C++ language.  In general, the accepted standard is **C++11**, which is supported by all major OSes and compilers currently in use.  Any future exceptions to this will be listed here, as well as in the "Compilation notes" section of the Colvars doc.

Because Colvars relies on the build system of each engine, there is nothing for you to do, whether you are using a precompiled build or making your own.  The only exception to this is VMD, whose build system requires some patching in order to work on modern architectures (details below).


### Backward compatibility notes

If you are reading this page following an error message from Colvars, you are probably using an earlier version of the code, which at the time was supporting C++11 as an option but not *requiring* it yet.  Consider upgrading your engine and/or Colvars, as you will also get new features and improvements.

We no longer support using extremely old compilers and OSes that are well past their end of life (for example, GCC 4.4 and CentOS 6).


### Status of C++ language support in MD engines (as of 2023-06-22)

- **GROMACS** follows the C++17 standard.

- **LAMMPS** follows the C++11 standard.

- **NAMD** follows the C++11 standard.

- **VMD** precompiled builds *do not* use C++11 yet, and the optional `CXX11` flag in VMD's `configure` script does not work yet.  Please apply the `c++11` patch from [this repository](https://github.com/giacomofiorin/vmd-patches/) to allow building VMD with C++11 whenever supported.


### Building VMD is quite complicated: could you just provide a working VMD build?

In theory we *could*, but due to UIUC's restrictive license we are curently unable to distribute any version of VMD, modified or not.  We hope that this restriction is lifted in the future, through the adoption of best practices for [research software sharing](https://datascience.nih.gov/tools-and-analytics/best-practices-for-sharing-research-software-faq).
