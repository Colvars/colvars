# Compute the reduced mass of a collective variable
# and write its trajectory to a file

# M = (\nabla xi m^-1 \nabla_xi)^-1
# where m is the diagonal matrix of Cartesian coordinate (ie atom) masses

proc cv_mass_traj { cv filename } {

  # Get atom ids and masses
  cv colvar $cv set collect_atom_ids 1
  set ids [cv colvar $cv getatomids]
  set sel [atomselect top "index $ids"]
  set inv_m [list]
  foreach mi [$sel get mass] {
    set inv [expr 1.0 / $mi]
    # one inverse mass per coordinate
    lappend inv_m $inv $inv $inv
  }

  cv colvar $cv set collect_gradient 1

  set out [open $filename "w"]
  set nf [molinfo top get numframes]

  for { set f 0 } { $f < $nf } { incr f } {
    if {[expr {$f % ($nf/20)}] == 0} { puts stdout "frame $f"}
    cv frame $f
    cv colvar $cv update
    set g [join [cv colvar $cv getgradients]]
    set M [expr 1.0 / [vecsum [vecmul [vecmul $g $inv_m] $g]]]
    puts $out "$f $M"
  }

  close $out
  puts "Data written to $filename"
  return
}
