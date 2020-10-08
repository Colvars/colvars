# -*- sh-basic-offset: 2; sh-indentation: 2; -*-

# Generated with Gromacs version 2020.3

# Command line to create the tpr files.
gmx_d grompp -f system.mdp -c system.gro -p system.top -o test.tpr

# Restart tpr
gmx_d convert-tpr -s test.tpr -o test.restart.tpr -nsteps 100

