# Write a VMD selection into a GROMACS index file.

# 1st argument is either a Tcl file channel or a file name: in the latter
# case, content will be appended to that file.
# 2nd argument is an atom selection proc, as returned by the atomselect
# command.
# 3rd argument is the name of the group.

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
