set alpha [atomselect top "alpha"]
set random_coords [list]
for { set i 0 } { ${i} < [${alpha} num] } { incr i } {
    set random_x [list]
    for { set d 0 } { ${d} < 3 } { incr d } {
        set r1 [expr rand()]
        set r2 [expr rand()]
        lappend random_x [expr ${r1} * (${r2} - 0.5)/abs(${r2} - 0.5)]
    }
    lappend random_coords ${random_x}
}

${alpha} set [list "x" "y" "z"] ${random_coords}
puts [measure center ${alpha}]
${alpha} writexyz random_vec_10.xyz
${alpha} delete
