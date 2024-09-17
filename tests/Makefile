.PHONY: all functional gromacs lammps namd


all: functional lammps namd # gromacs


functional:
	./build_tests.sh input_files/

# GROMACS tests need manual adjustments
# gromacs:
# 	./build_tests.sh ../gromacs/tests/library/

lammps:
	./build_tests.sh ../lammps/tests/library/

namd:
	./build_tests.sh ../namd/tests/library/

