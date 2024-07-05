# Create colvar for native contacts for current structure
# needs topology info, ie. resids
# assuming that resids are unique

proc native_contacts_cfg { {molid "top"} {seltext "alpha"} {cutoff 5} } {

    # no contacts up to third neighbor residue
    set min_res_dist 3

    set sel [atomselect $molid $seltext]
    set contacts [measure contacts $cutoff $sel]

    # Build lookup table of resids
    set res [dict create]
    foreach a [$sel list] r [$sel get resid] {
        dict set res $a $r
    }

    set nonlocal_a1 [list]
    set nonlocal_a2 [list]
    # contacts is a list of two lists of atom indices
    foreach a1 [lindex $contacts 0] a2 [lindex $contacts 1] {
        set seq_dist [expr abs([dict get $res $a2] - [dict get $res $a1])]
        if { $seq_dist > $min_res_dist } {
            lappend nonlocal_a1 $a1
            lappend nonlocal_a2 $a2
        }
    }

    set n_pairs [llength $nonlocal_a1]
    if { $n_pairs == 0} { return "" }
    set coeff [expr {1. / $n_pairs}]

    set cfg "colvar {
    name native_contacts
"
    foreach a1 $nonlocal_a1 a2 $nonlocal_a2 {
        # donor and acceptor are actually symmetric
        append cfg "
    hBond {
        componentCoeff $coeff
        donor $a1
        acceptor $a2
        cutoff $cutoff
    }"
    }
    append cfg "\n}\n"
    return $cfg
}

cv colvar native_contacts delete
cv config [native_contacts_cfg 0 "name CA" 7]
