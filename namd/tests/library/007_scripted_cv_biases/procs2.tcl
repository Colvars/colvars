
proc calc_Id {x} {
  return $x
}

proc calc_Id_gradient {x} {
  return 1.
}

proc calc_vector { x1 x2 } {
  return [concat [veclength $x2] $x1]
}

proc calc_vector_gradient { x1 x2 } {
  set c [expr 1. / [veclength $x2]]
  foreach { ex ey ez } [vecscale $x2 $c] {}
  return [list "\
  0. 0. 0. 0.\
  1. 0. 0. 0.\
  0. 1. 0. 0.\
  0. 0. 1. 0.\
  0. 0. 0. 1."\
 "$ex $ey $ez\
  0. 0. 0.\
  0. 0. 0.\
  0. 0. 0.\
  0. 0. 0."]
}

proc calc_colvar_forces { ts } {
  # Scripted harmonic bias on first component of colvar vec
  set center 18.
  set k 500.

  set d [lindex [cv colvar vec value] 0]
  set f [expr {$k * ($center - $d)}]
  cv colvar vec addforce [list $f 0. 0. 0. 0.]

  # Constant force on orientation part
  cv colvar vec addforce "0. 0. 100000. 0. 0."
  return
}
