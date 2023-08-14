#!/bin/csh -e
 
setenv gmx ~/local/gromacs-2021.6/bin/gmx
#setenv gmx ~/local/gromacs-2021-6-debug/bin/gmx

#setenv flag '-ntomp 1 -dlb yes'
setenv flag '-ntomp 1' ; # -ntmpi 1'

setenv bias_flag '-colvars colvars_ae_abf.dat'
#setenv bias_flag '-colvars colvars_nn_dihedrals_abf.dat'
#setenv bias_flag '-colvars colvars_abf.dat'

${gmx} mdrun $flag -deffnm md ${bias_flag}
