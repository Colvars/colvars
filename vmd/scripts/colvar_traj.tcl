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
