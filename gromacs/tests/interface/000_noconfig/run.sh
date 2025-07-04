# -*- sh-basic-offset: 2; sh-indentation: 2; -*-
TMPDIR=/tmp

# if script called from `run_tests.sh` in the parent folder, retrieve the correct command binary from the argument
if [ -x "$1" ]
then
    BINARY=$1
else
    BINARY="gmx_d"
fi

process_output () {
    output=${basename}.part0002
    for file in ${output}.* ; do
        # Skip if no files match (avoid literal '*.')
        [ -e "$file" ] || continue
        # Remove the part number
        mv -f ${file} ${file/.part0002/}
    done

    # Filter out the version numbers to allow comparisons
    grep "^colvars:" ${basename}.log \
    | grep -v 'Initializing the collective variables module' \
    | grep -v 'Using GROMACS interface, version' > ${basename}.colvars.out
    if [ -f ${basename}.colvars.state ] ; then
    grep -sv 'version' ${basename}.colvars.state \
        > ${TMPDIR}/${basename}.colvars.state.stripped && \
        mv -f ${TMPDIR}/${basename}.colvars.state.stripped ${basename}.colvars.state.stripped
    fi
}

options="-ntomp 2"

basename=test
$BINARY grompp -f nocvconfig.mdp -c ../Common/system.gro -p ../Common/system.top -o ${basename}.tpr  2> ${basename}.grompp.err 1> ${basename}.grompp.out
$BINARY mdrun -s ${basename}.tpr -deffnm ${basename} ${options} 2> ${basename}.err 1> ${basename}.out
process_output

basename=test.restart
# Restart tpr
$BINARY convert-tpr -s ${basename%.restart}.tpr -nsteps 100 -o ${basename}.tpr 2> ${basename}.grompp.err 1> ${basename}.grompp.out
$BINARY mdrun -s ${basename}.tpr -deffnm ${basename} -noappend -cpi ${basename%.restart}.cpt ${options} 2> ${basename}.err 1> ${basename}.out
process_output


rm -f *.tpr
