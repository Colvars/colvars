#!/bin/csh -e
 
#SBATCH --job-name=gro 				# replace the name
#SBATCH --mem=2048                              # replace with amount suitable for your job
#SBATCH --time=00-02:00:00                      # replace with amount suitable for your job
#SBATCH --partition=gpu
#SBATCH -w htc-gpu017
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --gres=gpu:2

setenv gmx gmx_mpi

#setenv flag '-ntomp 1 -dlb yes'
setenv flag '-ntomp 1' ; # -ntmpi 1'
setenv bias_flag ; #'-colvars colvars-abf.dat'

# Prepare run
${gmx} grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr

mpirun -np 1 ${gmx} mdrun $flag -deffnm md ${bias_flag}
