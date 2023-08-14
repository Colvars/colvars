#!/bin/csh -e
 
#SBATCH --job-name=gro 				# replace the name
#SBATCH --mem=2048                              # replace with amount suitable for your job
#SBATCH --time=00-02:00:00                      # replace with amount suitable for your job
#SBATCH --partition=gpu_micro
#SBATCH -w htc-gpu017
#SBATCH --nodes=1 
#SBATCH --gres=gpu:1

setenv pdbfile ../../pdb/ad.pdb

setenv flag '-ntomp 1 -ntmpi 1'

setenv gmx gmx


echo "Build system:"

echo "6" | $gmx pdb2gmx -f $pdbfile -water tip3p -o top.gro 

echo "create box:"
# Create box
${gmx} editconf -f top.gro -o box_gro -c -d 1.2 -bt cubic 

# Add Water
echo "Add water:" 
${gmx} solvate -cp box_gro -cs spc216.gro -o solv.gro -p topol.top 

# ION
echo "Add ion:" 
${gmx} grompp -f minim.mdp -maxwarn 4 -c solv.gro -p topol.top -o ions.tpr 
echo "13" | ${gmx} genion -s ions.tpr -o ions.gro -p topol.top -pname NA -nname CL -neutral 

# Energy minimization
${gmx} grompp -f minim.mdp -c ions.gro -p topol.top -o em.tpr
${gmx} mdrun $flag -deffnm em

# NVT 
${gmx} grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
${gmx} mdrun $flag -deffnm nvt

# NPT 
${gmx} grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
${gmx} mdrun $flag -deffnm npt

