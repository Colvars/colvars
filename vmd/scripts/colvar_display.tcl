###########################################################
# Real-time display of colvar values in a VMD scene
# with a graphical representation for rotation quaternions
# can be updated by pressing 'u' if coordinates change
#
# Example usage:
#
# cv molid top
# cv configfile my_colvars.in
#
# start_colvar_display "cv1 cv2" yellow   # text color
# set_colvar_display position { 0 0 5 }   # change text position
# stop_colvar_display
#
# start_rotation_display orientation
# stop_rotation_display
#
###########################################################

proc start_colvar_display { {cvs all} {text_color blue} {text_position {10 0 10}} {molid top} } {

  global text_pos text_col cv_handle cv_list vmd_frame

  set text_pos $text_position
  set text_col $text_color
  set cv_handle -1
  set cv_list $cvs

  if {! [string compare $molid top]} {
    set molid [molinfo top]
  }

  # draw values right away
  update_colvar_display vmd_frame $molid w

  # set trace
  trace add variable vmd_frame($molid) write update_colvar_display
  # set hotkey "u" to update the display
  user add key u "update_colvar_display vmd_frame $molid w"

  puts "Press 'u' to update the displayed values after modifying coordinates"
  puts "Remove the display with 'stop_colvar_display'"
  return
}


###########################################################
proc stop_colvar_display { {molid top} } {
  global cv_handle vmd_frame

  if { $cv_handle > -1 } {
    graphics $molid delete $cv_handle
  }
  if {! [string compare $molid top]} {
    set molid [molinfo top]
  }
  trace remove variable vmd_frame($molid) write update_colvar_display
  return
}

###########################################################
proc set_colvar_display { property value } {
  global text_pos text_col

  switch  $property  {
    "color"	{ set text_col $value }
    "position"	{ set text_pos $value }
    default	{ puts "Possible properties: color, position" }
  }

  update_colvar_display vmd_frame [molinfo top] w
  return
}


###########################################################
proc update_colvar_display {name molid op} {
  # name == "vmd_frame"
  # molid == molecule id of the newly changed frame
  # op == "w"

  global text_pos text_col cv_handle cv_list
  set f [molinfo $molid get frame]
  cv frame $f

  if {! [string compare $cv_list all]} {
    set cv [cv list]
  } else {
    set cv $cv_list
  }

  if { $cv_handle > -1 } {
    graphics $molid replace $cv_handle
  }

  foreach v $cv {
    # Calculate only requested variables
    set value [cv colvar $v update]
    # Parse vectors
    set value [regsub -all {,|\(|\)} $value ""]
    append t [format "%s=" $v]
    foreach x $value { append t [format "%.2f " $x] }
  }

  graphics $molid color $text_col
  set cv_handle [graphics $molid text $text_pos $t size 0.9 thickness 2]
  return
}

###########################################################
###########################################################
#   Display rotations as axis + "clock hands"
###########################################################
###########################################################
proc start_rotation_display { cv {color yellow} {position {0 0 0}} {molid top} } {

  global rot_pos rot_col rot_handle rot_cv vmd_frame

  set rot_pos $position
  set rot_col $color
  set rot_handle ""
  set rot_cv $cv

  if {! [string compare $molid top]} {
    set molid [molinfo top]
  }

  # draw values right away
  update_rotation_display vmd_frame $molid w

  # set trace
  global vmd_frame
  trace add variable vmd_frame($molid) write update_rotation_display
  return
}

###########################################################
proc stop_rotation_display { {molid top} } {
  global rot_handle vmd_frame

  foreach h $rot_handle {
    graphics $molid delete $h
  }
  if {! [string compare $molid top]} {
    set molid [molinfo top]
  }
  trace remove variable vmd_frame($molid) write update_rotation_display
  return
}

###########################################################
proc set_rotation_display { property value } {
  global rot_pos rot_col

  switch  $property  {
    "color"	{ set rot_col $value }
    "position"	{ set rot_pos $value }
    default	{ puts "Possible properties: color, position" }
  }

  update_rotation_display vmd_frame [molinfo top] w
  return
}


###########################################################
proc update_rotation_display {name molid op} {
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


  global rot_pos rot_col rot_handle rot_cv M_PI
  set f [molinfo $molid get frame]
  cv frame $f

  foreach h $rot_handle {
    graphics $molid delete $h
  }

  set value [cv colvar $rot_cv update]
  # Parse vector
  set value [regsub -all {,|\(|\)} $value ""]
  foreach {q0 q1 q2 q3} $value {}

  graphics $molid color $rot_col

  set theta [expr {2.0 * 180. / $M_PI * acos($q0)}]
  if {$theta*$theta < 1e-12} {
      lappend  rot_handle [graphics $molid sphere $rot_pos radius $R resolution $res]
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
      lappend  rot_handle [graphics $molid cylinder [vecsub $rot_pos $u] [vecadd $rot_pos $u] radius $R resolution $res filled yes]
      set tip [vecadd [vecadd $rot_pos $u] [vecscale $u [expr {2. * $tip_ratio * $R / $L}]]]
      lappend  rot_handle [graphics $molid cone [vecadd $rot_pos $u] $tip radius [expr {1.5 * $R}] resolution $res]
      lappend  rot_handle [graphics $molid cylinder $rot_pos $e1 radius $r resolution $res filled yes]
      lappend  rot_handle [graphics $molid cylinder $rot_pos $e2 radius $r resolution $res filled yes]
  }

  return
}

