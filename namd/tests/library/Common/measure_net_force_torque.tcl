# Insert into a NAMD script to get output with the total applied force and torque
# -*- mode: tcl; tcl-indent-level: 2; -*-

if { [info procs num_atoms] == "" } {
    proc num_atoms {} {
        return 104
    }
}

tclForces on
tclForcesScript {
  set t 0
  for {set i 1} {$i <= [num_atoms]} {incr i} {
    addatom $i
  }
  proc calcforces {} {
    global t
    global prev_c

    if { $t == 0} {
      loadcoords prev_c
    } else {
      loadforces f
      loadcoords c
      set tot [list 0. 0. 0.]
      set torque [list 0. 0. 0.]

      for {set i 1} {$i <= [num_atoms]} {incr i} {
        if { [info exists f($i)] } {
          set tot [vecadd $tot $f($i)]

          foreach {x y z} $prev_c($i) {}
          foreach {fx fy fz} $f($i) {}
          set to [list [expr $y*$fz-$z*$fy] [expr $z*$fx-$x*$fz] [expr $x*$fy-$y*$fx]]
          set torque [vecadd $torque $to]
        }
      }
      #      print "[expr $t-1]   Force: ( $tot )   Torque: ( $torque )"
      print "TOTAL FORCE: ( $tot )"
      print "TOTAL TORQUE: ( $torque )"
      array set prev_c [array get c]
    }
    incr t
    return
  }
}

