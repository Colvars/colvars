
proc calc_Id { x } {
  return $x
}
proc calc_Id_gradient { x } {
  return 1.
}
proc calc_Id_quaternion { x } {
  return $x
}
proc calc_Id_quaternion_gradient { x } {
  return [list "(1., 1., 1., 1.)"]
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
