
# Read a GROMACS-style .ndx file, associate its groups to VMD atom selection macros
# e.g. the group "Protein" can be used as: atomselect top "ndx_Protein"
# See colvartools/write_index_group.tcl for the reverse

proc read_index_file { filename { prefix "ndx_" } { suffix "" } } {
    array set index_groups {}
    set input [open ${filename} "r"]
    set group_name ""
    set line_number 0
    while { [gets ${input} line] >= 0 } {
        set words [list {*}${line}]
        incr line_number
        if { [llength ${words}] == 0 } {
            continue
        }
        if { [lindex ${words} 0] == "\[" } {
            if { [lindex ${words} 2] != "\]" } {
                puts stderr "ERROR: Wrong group header at line ${line_number} in ${filename}"
                ${input} close
                return 1
            }
            set group_name [lindex ${words} 1]
            set index_groups(${group_name}) [list]
        } else {
            lappend index_groups(${group_name}) {*}${words}
        }
    }
    foreach group_name [lsort [array names index_groups]] {
        set indices [lindex [array get index_groups ${group_name}] 1]
        if { [llength ${indices}] > 0 } {
            puts "Index group \"${group_name}\": [llength ${indices}] atoms -> atomselect macro ${prefix}${group_name}${suffix}"
            atomselect macro ${prefix}${group_name}${suffix} "(serial [join ${indices}])"
        }
    }
    array unset index_groups
    close ${input}
}
