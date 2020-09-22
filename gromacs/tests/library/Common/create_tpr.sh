# -*- sh-basic-offset: 2; sh-indentation: 2; -*-

# Generated with Gromacs version 2020.3

# Note: The peptide coordinates have to be changed compared to namd.
# We need to center the peptides inside the box (editconf) to avoid some arbitrary translations made
# by gromacs mdrun at the first step, due to the fact the peptide is at a corner of the box (without centering).
# This change a bit (at 1e-5 precision) the results of some tests.

# Command line to create the topologies files from da.coor.pdb
# Termini  are NH2 and CT2
# Forcefield is Charmm22 (Charmm27 without CMAP)
echo "1 2"|gmx_d pdb2gmx -f da.coor.pdb -ff charmm27 -water none -ter -nocmap -ignh -o da_nobox.pdb -p da.top

# !! MANUALLY CHANGES !!
# 1. In da_nobox.pdb, replace H coordinates with the ones in da.coor.pdb
# 2. In da.top, change the charges of the 5 first atoms of ALA1 due to a bad psf in the namd tests
# New charges in the order : -0.62 0.31 0.31 -0.10 0.10
# 3. Change also the type of the 3 first atoms. New types in order: NH3,HC,HC

# Generate box
# Center the peptide inside the box (cf Note)
gmx_d editconf -f da_nobox.pdb -bt cubic -d 1.5 -o da.pdb

# Generate trr for starting coordinates and velocities
# Pick up the box dimension from editconf
python create_trr.py

# Command line to create the tpr files.
gmx_d grompp -f test.mdp -c da.pdb -p da.top -t da.trr -o test
# Restart tpr
gmx_d convert-tpr -s test.tpr -nsteps 40 -o test.restart.tpr

