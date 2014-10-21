
proc calc_Id { x } {
  return $x
}
proc calc_Id_gradient { x } {
  return 1.
}
proc calc_Id_pair { x1 x2 } {
  return $x2
}
proc calc_Id_pair_gradient { x1 x2 } {
  return [list {(0., 0., 0., 0.)} 1.]
}

proc calc_colvar_forces { ts } {
  # Scripted harmonic bias on colvar ds
  set center 18.
  set k 500.

  set ds [colvars colvar ds value]
  set f [expr {$k * ($center - $ds)}]
  set r [colvars colvar ds addforce $f]
  return
}
