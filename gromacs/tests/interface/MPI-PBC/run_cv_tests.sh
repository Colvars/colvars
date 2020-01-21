
# Disabling GPU for reproducibility

gmx_d mdrun -s test.tpr -nb cpu  -colvars cv_nobias.colvars.dat
mkdir cv_nobias
rm -f cv_nobias/*
mv md.* cv_nobias/
mv *.xtc cv_nobias/

mpirun -np 4 gmx_mpi_d mdrun -s test.tpr -nb cpu  -colvars cv_nobias.colvars.dat
mkdir cv_nobias_mpi
rm -f cv_nobias_mpi/*
mv md.* cv_nobias_mpi/
mv *.xtc cv_nobias_mpi/

gmx_d mdrun -s test.tpr -nb cpu  -colvars cv_bias.colvars.dat
mkdir cv
rm -f cv/*
mv md.* cv/
mv *.xtc cv/

mpirun -np 4 gmx_mpi_d mdrun -s test.tpr -nb cpu  -colvars cv_bias.colvars.dat
mkdir cv_mpi
rm -f cv_mpi/*
mv md.* cv_mpi/
mv *.xtc cv_mpi/

# Compare MPI and non-MPI simulations
spiff cv_nobias/md.colvars.traj cv_nobias_mpi/md.colvars.traj 
spiff cv/md.colvars.traj  cv_mpi/md.colvars.traj
