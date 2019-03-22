# colvar_traj: prints the complete colvar trajectory for the loaded
# trajectory
# save_colvar_traj_file <file>: saves the trajectory to given file name


proc colvar_traj {} {

  set s [cv printframelabels]

  set nf [molinfo top get numframes]
  for {set f 0} {$f< $nf} {incr f} {
    cv frame $f
    cv update
    append s [cv printframe]
  }
  return $s
}

proc colvar_save_traj_file { {fileName vmd.colvars.traj} } {

  puts "Writing colvars trajectory to file $fileName"
  set o [open $fileName w]
  puts $o [colvar_traj]
  close $o
}

proc colvar_plot_traj {} {

  set nf [molinfo top get numframes]
  set x [list]
  for {set f 0} {$f< $nf} {incr f} { lappend x $f }

  foreach c [cv list] {
    set y($c) [list]
  }
  for {set f 0} {$f< $nf} {incr f} {
    cv frame $f
    cv update
    foreach c [cv list] {
      lappend y($c) [cv colvar $c value]
    }
  }
  set plothandle [multiplot -title "Colvars trajectory" -xlabel "Frame" -ylabel "Colvar value"]
  foreach c [cv list] {
    $plothandle add $x $y($c) -legend $c
  }
  $plothandle replot
}
