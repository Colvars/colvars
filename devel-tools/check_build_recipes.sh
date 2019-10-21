#!/bin/bash

for f in src/*.cpp ; do
    f=$(basename $f)
    namd_target="\$(DSTDIR)/${f%.cpp}.o"
    if ! grep -q ${namd_target} namd/colvars/src/Makefile.namd ; then
        echo "${namd_target} missing from namd/colvars/src/Makefile.namd"
    fi
    vmd_source=${f%.cpp}.C
    if ! grep -q ${vmd_source} vmd/src/colvars_files.pl ; then
        echo "${vmd_source} missing from vmd/src/colvars_files.pl"
    fi
    if ! grep -q ${f} lammps/lib/colvars/Makefile.common ; then
        echo "${f} missing from lammps/lib/colvars/Makefile.common"
    fi
done

