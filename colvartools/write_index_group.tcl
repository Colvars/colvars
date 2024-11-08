# Write a VMD selection into a GROMACS index file.

# Parameters:
# 1- either a Tcl file channel or a file name: in the latter
#   case, content will be appended to that file.
# 2- an atom selection proc, as returned by the atomselect
# command.
# 3- the name of the group.

proc write_index_group { ndxfile sel name } {
    # Check that the name does not contain spaces or tabs
    if { ([string first " " "${name}"] >= 0) ||
         ([string first "	" "${name}"] >= 0) } {
        puts stderr "Error: group name may not contain spaces or tabs."
        return
    }

    if { [catch [list tell "${ndxfile}"]] } {
        # ${ndxfile} is the filename
        set output [open "${ndxfile}" "a"]
    } else {
        # ${ndxfile} is the channel
        set output ${ndxfile}
    }

    puts ${output} [format "\[ %s \]" ${name}]
    set line_buf 0
    set atom_count 0
    foreach num [${sel} get serial] {
        incr atom_count
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

    if { "${output}" != "${ndxfile}" } {
        close ${output}
    }
}

# Write a GROMACS index file suitable for computing the alpha-helix
# content of a helical segment

# Parameters:
# 1- either a Tcl file channel or a file name: in the latter
#   case, content will be appended to that file.
# 2- a selection text returning contiguous amino-acid residues
# 3- optional: molecule id (default: top)
# 4- optional: prefix for the group names (default: alpha_)

proc write_alpha_groups { ndxfile seltext { mol top } { prefix alpha_ } } {

    foreach atomname { N CA O } {
        set sel [atomselect $mol "($seltext) and name $atomname"]
        write_index_group $ndxfile $sel "${prefix}${atomname}"
        $sel delete
    }
}

# Write a GROMACS index file suitable for computing the dihedralPC
# projection of a peptide chain

# Parameters:
# 1- either a Tcl file channel or a file name: in the latter
#   case, content will be appended to that file.
# 2- a selection text returning contiguous amino-acid residues
# 3- optional: molecule id (default: top)
# 4- optional: prefix for the group names (default: alpha_)

proc write_dihedralPC_groups { ndxfile seltext { mol top } { prefix dihed_ } } {

    foreach atomname { CA N C } {
        set sel [atomselect $mol "($seltext) and name $atomname"]
        write_index_group $ndxfile $sel "${prefix}${atomname}"
        $sel delete
    }
}