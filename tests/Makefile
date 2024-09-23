.PHONY: all functional gromacs lammps namd


all: functional gromacs lammps namd

functional:
	./build_tests.sh input_files/

gromacs:
	./build_tests.sh ../gromacs/tests/library/

lammps:
	./build_tests.sh ../lammps/tests/library/

namd:
	./build_tests.sh ../namd/tests/library/

