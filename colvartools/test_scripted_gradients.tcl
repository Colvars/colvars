# Numerical test for analytical scripted gradients
#
# example usage in tclsh:
#
# source test_scripted_gradients.tcl
# source ../namd/tests/library/007_scripted_cv_biases/procs2.tcl 
# test_grad vector {{0. 1. 0 0} {1. 2. 3. 4. 5. 6.}}

if { [info tclversion] < 8.5 } {
    puts "###### Sorry, this script uses the {*} idiom of Tcl 8.5 and later. Please update. ######"
}

# Utility functions

# MSD between two lists of lists
proc msd {vvx vvy} {
    set s 0.
    foreach vx $vvx vy $vvy {
        foreach x $vx y $vy {
            set s [expr {$s + ($x-$y)*($x-$y)}]
        }
    }
    return $s
}

proc vecsub {a b} {
    set r {}
    foreach x $a y $b {lappend r [expr {$x - $y}]}
    return $r
}
proc veclength2 {v} {
    set l2 0.
    foreach x $v { set l2 [expr {$l2 + $x*$x}]}
    return $l2
}
proc veclength {v} {
    set l2 0.
    foreach x $v { set l2 [expr {$l2 + $x*$x}]}
    return [expr {sqrt($l2)}]
}
proc vveclength2 {vv} {
    set l2 0.
    foreach v $vv {
       foreach x $v {set l2 [expr {$l2 + $x*$x}]}
    }
    return $l2
}
proc vecscale {v s} {
    if {[llength $s]>1} {
        set t $v
        set v $s
        set s $t
    }
    set r {}
    foreach x $v { lappend r [expr {$s * $x}]}
    return $r
}

# Gradient from finite differences

proc grad_fdiff {func x0 dx {widths {}}} {

    # TODO normalize dx by widths
    set ncvc [llength $x0]
    set df [llength [calc_${func} {*}$x0]]

    set g_fdiff {}
    # gradient wrt each cvc
    for {set i 0} {$i < $ncvc} {incr i} {
        set xi0 [lindex $x0 $i]
        set d [llength $xi0]
        for {set k 0} {$k<$df} {incr k} {set gi($k) {}}

        # gradient wrt each cvc's components
        for {set j 0} {$j < $d} {incr j} {
            set x1 [lreplace $x0 $i $i [lreplace $xi0 $j $j [expr [lindex $xi0 $j] - $dx/2.]]] 
            set f1 [calc_${func} {*}$x1]
            set x2 [lreplace $x0 $i $i [lreplace $xi0 $j $j [expr [lindex $xi0 $j] + $dx/2.]]] 
            set f2 [calc_${func} {*}$x2]

            # the function value is a vector
            set k 0
            foreach a $f1 b $f2 {
                lappend gi($k) [expr {($b - $a) / $dx}]
                incr k
            }
        }
        # Need to transpose the result
        set tmp {}
        for {set k 0} {$k<$df} {incr k} {
            set tmp [concat $tmp $gi($k)]
        }
        lappend g_fdiff $tmp
    }
    return $g_fdiff
}

# User front-end
# Outputs dx, relative RMSD for given dx
# for decreasing values of dx

proc test_grad {func x0 {dx 1.} {widths {}}} {
    # call the value function first: might be necessary to initialize internal variables
    # for the gradient calculation
    calc_${func} {*}$x0

    set ga [calc_${func}_gradient {*}$x0]
    for {set n 0} {$n < 5} {incr n} {
        set gd [grad_fdiff $func $x0 $dx $widths]
        puts "$dx [expr sqrt([msd $ga $gd] / [vveclength2 $ga])]"
        set dx [expr {$dx / 10.}]
    }
}
