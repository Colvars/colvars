
# Regenerate only if necessary
# beware of random seed - gen-seed should be in mdp
#gmx_d grompp -f system.mdp  -p system.top   -o test.tpr  -c system.gro

# Disabling GPU (for nonbonded only?)

gmx_d mdrun  -s test.tpr -nb cpu
mkdir nocv
rm -f nocv/*
mv md.* nocv/
mv *.xtc nocv/

mpirun -np 4 gmx_mpi_d mdrun  -s test.tpr -nb cpu
mkdir nocv_mpi
rm -f nocv_mpi/*
mv md.* nocv_mpi/
mv *.xtc  nocv_mpi/


