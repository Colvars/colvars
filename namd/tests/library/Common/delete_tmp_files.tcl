if { [info exists mol_name] == 0 } {
    set mol_name "da"
}

if { ${mol_name} == "da" } {
    file delete -force "index.ndx"
    file delete -force "rmsd_atoms_refpos.xyz"
    file delete -force "rmsd_atoms_random.xyz"
    file delete -force "heavy_atoms_refpos.xyz"
}
