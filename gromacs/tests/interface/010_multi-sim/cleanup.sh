#!/bin/bash
# -*- sh-basic-offset: 2; sh-indentation: 2; -*-

# Script to remove simulations files
# Called by a `run_tests.sh` scripts

for i in a b c d
do
    cd $i
    rm -f *.xtc *.trr *.edr *.cpt *.gro *.log *~ \#*\# mdout.mdp
    rm -f *.state *.old *.traj
    cd ..
done
