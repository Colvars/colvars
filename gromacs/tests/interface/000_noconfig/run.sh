# -*- sh-basic-offset: 2; sh-indentation: 2; -*-
TMPDIR=/tmp

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

# For GROMACS 2024 with Colvars MdModule

options="-ntmpi 1 -ntomp 2"

basename=test
gmx_d grompp -f nocvconfig.mdp -c ../Common/system.gro -p ../Common/system.top -o ${basename}.tpr  2> ${basename}.grompp.err 1> ${basename}.grompp.out
gmx_d mdrun -s ${basename}.tpr -deffnm ${basename} ${options} 2> ${basename}.err 1> ${basename}.out
process_output

basename=test.restart
# Restart tpr
gmx_d convert-tpr -s ${basename%.restart}.tpr -nsteps 100 -o ${basename}.tpr 2> ${basename}.grompp.err 1> ${basename}.grompp.out
gmx_d mdrun -s ${basename}.tpr -deffnm ${basename} -noappend -cpi ${basename%.restart}.cpt ${options} 2> ${basename}.err 1> ${basename}.out
process_output


rm -f *.tpr
