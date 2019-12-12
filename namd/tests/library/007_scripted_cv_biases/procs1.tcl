
proc calc_colvar_forces { ts } {
  # Scripted constant force on orientation quaternion

  cv colvar o addforce "0. 100000. 0. 0."
  return
}
