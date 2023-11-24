# -*- sh-basic-offset: 2; sh-indentation: 2; -*-

# For GROMACS 2024 with Colvars MdModule

gmx_d grompp -f nocvconfig.mdp -c ../Common/system.gro -p ../Common/system.top -o test.tpr
# Restart tpr
gmx_d convert-tpr -s test.tpr -o test.restart.tpr -nsteps 100


basename=test
options="-ntmpi 1 -ntomp 2"

gmx_d mdrun -s test.tpr -deffnm ${basename} ${options} &> ${basename}.out

gmx_d mdrun -s test.tpr -deffnm ${basename} -noappend -cpi ${basename%.restart}.cpt -colvars_restart test.colvars.state.dat ${options} &> ${basename}.out
