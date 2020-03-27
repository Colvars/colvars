mol new da.psf waitfor all

foreach pdbfile [list da.pdb da-rotated.pdb] {
    mol addfile ${pdbfile} waitfor all
    set_elements
    set sel [atomselect top "(all)"]
    ${sel} set radius 2.0
    volmap \
        density \
        ${sel} \
        -weight [${sel} get atomicnumber] \
        -minmax { { -10.0 -10.0 -10.0} {  10.0  10.0  10.0} } \
        -res 1.0 \
        -allframes \
        -combine avg \
        -checkpoint 0 \
        -o [file rootname ${pdbfile}].density.dx
    ${sel} delete
    animate delete all
}

mol delete top
mol new test.pdb
set sel [atomselect top "(all)"]
${sel} set radius 2.0
volmap \
    density \
    ${sel} \
    -weight [expr sqrt(8.0*acos(-1)*acos(-1)*acos(-1))] \
    -minmax { { -10.0 -10.0 -10.0} {  10.0  10.0  10.0} } \
    -res 1.0 \
    -allframes \
    -combine avg \
    -checkpoint 0 \
    -o gaussian.density.dx
${sel} delete

quit
