proc write_index_group { output sel label { n_atoms_max 0 } } {
    if { [${sel} num] == 0 } return
    puts stderr "Writing index group \"${label}\": [${sel} num] atoms"
    puts -nonewline ${output} [format "\[ %s \]\n" ${label}]
    set line_buf 0
    set atom_count 0
    foreach num [${sel} get serial] {
        incr atom_count
        if { (${n_atoms_max} > 0) && (${atom_count} > ${n_atoms_max}) } break
        puts -nonewline ${output} [format " %9d" ${num}]
        set line_buf [expr ${line_buf} + 10]
        if { ${line_buf} > 70 } {
            set line_buf 0
            puts -nonewline ${output} "\n"
        }
    }
    if { ${line_buf} > 0 } {
        puts -nonewline ${output} "\n"
    }
    puts -nonewline ${output} "\n"
}


if { [info exists mol_name] == 0 } {
    set mol_name [file rootname [molinfo top get name]]
}


if { ${mol_name} == "da" } {

set protein_sel "((protein) or (segid BH))" 

set index_file [open "da.ndx" "w"]

write_index_group ${index_file} \
    [atomselect top "${protein_sel}"] \
    "Protein"
write_index_group ${index_file} \
    [atomselect top "(${protein_sel} and (noh))"]\
    "Protein_noH"
write_index_group ${index_file} \
    [atomselect top "(${protein_sel} and (backbone))"] \
    "Protein_Backbone"
write_index_group ${index_file} \
    [atomselect top "(${protein_sel} and (alpha))"] \
    "Protein_C-alpha"
write_index_group ${index_file} \
    [atomselect top "(${protein_sel} and (alpha))"] \
    "RMSD_atoms"

write_index_group ${index_file} \
    [atomselect top "(${protein_sel} and (alpha) and resid 1 2)"] \
    "Protein_C-alpha_1_2"
write_index_group ${index_file} \
    [atomselect top "(${protein_sel} and (alpha) and resid 9 10)"] \
    "Protein_C-alpha_9_10"

for { set resid 1 } { ${resid} <= 10 } { incr resid } {
    write_index_group ${index_file} \
        [atomselect top "(${protein_sel} and (alpha) and resid ${resid})"] \
        "Protein_C-alpha_${resid}"

    write_index_group ${index_file} \
        [atomselect top "(${protein_sel} and (name N CA C O) and resid ${resid})"] \
        "group${resid}"
}

write_index_group ${index_file} \
    [atomselect top "(${protein_sel} and (not hydrogen))"] \
    "heavy_atoms"

close ${index_file}

}

quit

