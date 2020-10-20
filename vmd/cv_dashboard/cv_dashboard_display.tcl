###################################################################
# Forked from colvar_display.tcl for integration into the Dashboard
#
# Only rotation_display, colvar_display is not integrated yet
#
###################################################################


###########################################################
###########################################################
#   Display rotations as axis + "clock hands"
###########################################################
###########################################################
proc ::cv_dashboard::start_rotation_display { cv {color yellow} {position {0 0 0}} {molid top} } {

  global rot_pos rot_col rot_handles rot_cv

  if { [llength $cv] != 1 } {
    tk_messageBox -icon error -title "Colvars Dashboard Error"\
      -message "Select exactly 1 rotation colvar to enable rotation display.\n"
    return
  }

  if { [info exists rot_cv] } {
    stop_rotation_display
  }

  set rot_pos $position
  set rot_col $color
  set rot_handles {}
  set rot_cv $cv

  if {! [string compare $molid top]} {
    set molid [molinfo top]
  }

  # draw values right away
  update_rotation_display
}

###########################################################
proc ::cv_dashboard::stop_rotation_display {} {
  global rot_handles rot_cv
  set molid $::cv_dashboard::mol

  foreach h $rot_handles {
    graphics $molid delete $h
  }

  unset rot_cv
}

###########################################################
proc ::cv_dashboard::set_rotation_display { property value } {
  global rot_pos rot_col

  switch  $property  {
    "color"	{ set rot_col $value }
    "position"	{ set rot_pos $value }
    default	{ puts "Possible properties: color, position" }
  }

  update_rotation_display
}


###########################################################
proc ::cv_dashboard::update_rotation_display {} {
  # name == vmd_frame
  # molid == molecule id of the newly changed frame
  # op == w

  # Geometric parameters: lengths and radii of cylinders
  set L 12.     ;# length of main axis
  set R 0.7     ;# radius of main axis
  set l 5.      ;# length of "clock handles"
  set r 0.5     ;# radius of "clock handles"
  set tip_ratio 4. ;# cone length divided by main radius
  set res 16    ;# resolution of all objects
  ######################################################

  global rot_pos rot_col rot_handles rot_cv M_PI

  # Nothing to update
  if { ![info exists rot_cv] } { return }

  set molid $::cv_dashboard::mol

  foreach h $rot_handles {
    graphics $molid delete $h
  }

  # Forget variables that have been deleted (*after* deleting the graphics above)
  if { [lsearch $::cv_dashboard::cvs $rot_cv] == -1 } {
    unset rot_cv
    return
  }

  # CV was already updated by update_frame
  set value [cv colvar $rot_cv value]
  foreach {q0 q1 q2 q3} $value {}

  graphics $molid color $rot_col

  set theta [expr {2.0 * 180. / $M_PI * acos($q0)}]
  if {$theta*$theta < 1e-12} {
      lappend  rot_handles [graphics $molid sphere $rot_pos radius $R resolution $res]
  } else {

      set factor [expr {$L / sqrt(1-$q0*$q0) / 2.}]
      # u has the half-length of the main cylinder
      set u [vecscale $factor "$q1 $q2 $q3"]
      if { $q1 == 0. } {
          set v1 [vecnorm [veccross $u "0. 1. 0."]]
      } else {
          set v1 [vecnorm [veccross $u "1. 0. 0."]]
      }
      set v2 [vectrans [transabout $u $theta] $v1]

      set e1 [vecadd $rot_pos [vecscale $l $v1]]
      set e2 [vecadd $rot_pos [vecscale $l $v2]]
      lappend  rot_handles [graphics $molid cylinder [vecsub $rot_pos $u] [vecadd $rot_pos $u] radius $R resolution $res filled yes]
      set tip [vecadd [vecadd $rot_pos $u] [vecscale $u [expr {2. * $tip_ratio * $R / $L}]]]
      lappend  rot_handles [graphics $molid cone [vecadd $rot_pos $u] $tip radius [expr {1.5 * $R}] resolution $res]
      lappend  rot_handles [graphics $molid cylinder $rot_pos $e1 radius $r resolution $res filled yes]
      lappend  rot_handles [graphics $molid cylinder $rot_pos $e2 radius $r resolution $res filled yes]
  }

  return
}

