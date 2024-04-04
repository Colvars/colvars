# Load energy difference values from a NAMD FEP output file
# To be used as scripted colvar

proc load_fepout { filenames } {

    set dE_back 0 ;# Constant to signal backward time step in IDWS

    set dE [list]
    set i 0
    foreach n $filenames {
        puts "Opening $n ..."
        set fd [open $n "r"]
        while { [gets $fd line] >= 0 } {
            if { [string range $line 0 3] == "FepE" } {
                if { [lindex $line 0] == "FepEnergy:" } {
                    lappend dE [lindex $line 6]
                } else {
                    # IDWS step, do not mix with fwd dE values
                    lappend dE $dE_back
                }
                incr i
                if [expr {$i % 100000 == 0}] { puts "Reading: $i"}
            }
        }
        close $fd
    }

    set ne [llength $dE]
    puts "Read $ne dE values"
    set nf [molinfo top get numframes]

    if { $nf > 0 && [expr {$ne % $nf}] == 0 } {
        set skip 0
        set stride [expr {$ne / $nf}]
        puts "Assuming $stride deltaE samples per trajectory frame."
    } elseif { $nf > 1 && [expr {$ne % ($nf - 1)}] == 0 } {
        set skip 1
        set stride [expr {$ne / $nf}]
        puts "Assuming $stride deltaE samples per trajectory frame, by skipping first frame."
        # Unset value at first frame (no info)
        $::atom0_sel frame 0
        $::atom0_sel set user 0
    } else {
        puts "ERROR: can't match $ne deltaE values with $nf frames."
        return
    }

    for {set f $skip} {$f < [molinfo top get numframes]} {incr f} {
        $::atom0_sel frame $f
        $::atom0_sel set user [lindex $dE [expr $stride * $f]]
    }
}

set atom0_sel [atomselect top "index 0"]

proc calc_user0 { x } {
    $::atom0_sel frame [cv frame]
    return [$::atom0_sel get user]
}

cv molid top
cv config "
colvar {
    name alch_dE
    scriptedFunction user0
    distance {
        # The distance component is just a placeholder
        group1 {
            dummyatom (0, 0, 0)
        }
        group2 {
            dummyatom (0, 0, 0)
        }
    }
}
"