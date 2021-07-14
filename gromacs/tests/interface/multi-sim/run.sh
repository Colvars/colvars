# Attempt replica exchange every ten steps
mpirun  -np 16 gmx_mpi_d  mdrun -multidir a b c d -s test.tpr  -deffnm test -colvars ../test.dat -replex 10
