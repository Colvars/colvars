gmx_d grompp -f system.mdp  -p system.top   -o test.tpr  -c system.gro

gmx_d mdrun -s test.tpr -nb cpu  -colvars cv_bias.colvars.dat -deffnm md 

ln -s md.colvars.state md.colvars.state.dat
# Extend tpr
gmx_d convert-tpr -s test.tpr -o test_restart.tpr -nsteps 1000

gmx_d mdrun -s test_restart.tpr -nb cpu  -colvars cv_bias.colvars.dat -deffnm restart -cpi state.cpt -noappend -colvars_restart md.colvars.state.dat

