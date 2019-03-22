#############################################################
# Namespace for global parameters
# Set lambda, freq, and tolerance eg. in NAMD config file
# WARNING: coordinate z cannot be used without s, as it reuses
# the image flags set by s
# Furthermore, they have to be defined in that order: s then z
# See example file pathCV_example.namd
# Contributed by Christophe Chipot <chipot@ks.uiuc.edu>
# and Jérôme Hénin <henin@ibpc.fr>
#############################################################

namespace eval pathCV {
    set step      0
}

#############################################################
# s(X)
#############################################################

proc calc_pathCVs { args } {

    upvar pathCV::step step
    upvar pathCV::lambda lambda
    upvar pathCV::tolerance tolerance
    upvar pathCV::freq freq
    upvar pathCV::min_images min_images
    upvar pathCV::u u
    upvar pathCV::v v
    upvar pathCV::flags flags
    upvar pathCV::flags_next flags_next

    incr step

    if { [ expr {($step % $freq) == 0 }]} {
        # Enable all cvcs to allow for selection on next timestep
        set flags_next { }
        foreach x $args {
           lappend flags_next 1
        }
        cv colvar s cvcflags $flags_next
    }

    if { [ expr {($step % $freq) == 1}] } {
        # Flags are all set to one from previous timestep
        set flags {}
        foreach x $args {
           lappend flags 1
        }
        # Select relevant cvcs with new flags
        set flags_next {}
        set n_nonzero 0
        foreach x $args {
           set expo [expr {exp(- $lambda * $x * $x)}]
           if { $expo < $tolerance } {
              set flag 0
           } else {
              set flag 1
           }
           incr n_nonzero $flag
           lappend flags_next $flag
        }
        if { $n_nonzero < $min_images } {
          puts "PathCV WARNING: less than $min_images RMSDs are above threshold, selecting $min_images closest images"
          set sorted [lsort -real $args]
          set max [expr {[lindex $sorted [expr {$min_images - 1}]] + 1e-12}]
          set flags_next {}
          foreach x $args {
            if { $x < $max } {
              lappend flags_next 1
            } else {
              lappend flags_next 0
            }
          }
        }
        cv colvar s cvcflags $flags_next
        puts "STEP $step   next pathCV flags  $flags_next"
    }

    if { [ expr {$step % $freq} ] == 2 } {
        # Only use CVCs selected on previous step
        set flags $flags_next
    }

    set N [llength $args]

    set i 0
    set u 0.0
    set v 0.0
    foreach x $args c $flags {
        set u [expr {$u + $i * $c * exp(- $lambda * $x * $x)}]
        set v [expr {$v + $c * exp(-$lambda * $x * $x)}]
        incr i
    }

    set s [expr {1.0 / ($N - 1.0) * $u / $v}]

    return $s
}


#############################################################
# ds(X)/dX
#############################################################

proc calc_pathCVs_gradient { args } {

    upvar pathCV::lambda lambda
    upvar pathCV::u u
    upvar pathCV::v v
    upvar pathCV::flags flags

    set grad {}
    set N [llength $args]

    set i -1
    foreach x $args c $flags {
       incr i
       if { $c == 0 } { continue }
       set uprime [expr {-2.0 * $i * $lambda * $x * exp(- $lambda * $x * $x)}]
       set vprime [expr {-2.0 * $lambda * $x * exp(- $lambda * $x * $x)}]

       set ds [expr {1.0 / ($N - 1.0) * ($uprime * $v - $vprime * $u) / ($v * $v)}]
       lappend grad $ds
    }

    return $grad
}


#############################################################
# z(X)
#############################################################

proc calc_pathCVz { args } {

    upvar pathCV::lambda lambda
    upvar pathCV::step step
    upvar pathCV::tolerance tolerance
    upvar pathCV::freq freq
    upvar pathCV::sum sum
    upvar pathCV::flags flags
    upvar pathCV::flags_next flags_next

    set N [llength $args]
    if { [llength $flags] != $N } {
        puts "pathCV z ERROR: wrong number of CVC flags ([llength $flags]) for given RMSD values ($N)"
    }

    if { [ expr {$step % $freq <= 1}]} {
        # Re-using the flags calculated by colvar s
        # assuming that s is defined before z, and hence updated first
        cv colvar z cvcflags $flags_next
    }

    set sum 0.0
    foreach x $args c $flags {
        if { $c == 0 } { continue }
        set sum [expr {$sum + exp(- $lambda * $x * $x)}]
    }
    return [expr {- 1.0 / $lambda * log($sum)}]
}


#############################################################
# dz(X)/dX
#############################################################

proc calc_pathCVz_gradient { args } {

    upvar pathCV::lambda lambda
    upvar pathCV::step step
    upvar pathCV::tolerance tolerance
    upvar pathCV::freq freq
    upvar pathCV::sum sum
    upvar pathCV::flags flags

    set N [llength $args]
    if { [llength $flags] != $N } {
        puts "pathCV z gradient ERROR: wrong number of CVC flags ([llength $flags]) for given RMSD values ($N)"
    }

    set grad {}
    foreach x $args c $flags {
        if { $c == 0 } { continue }
        set sum_prime [expr {-2.0 * $lambda * $x * exp(- $lambda * $x * $x)}]
        set dz [expr -1.0 / $lambda * $sum_prime / $sum]
        lappend grad $dz
    }

  return $grad
}

