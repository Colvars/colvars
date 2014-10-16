# Real-time display of colvar values in a VMD scene
#
# Example usage:
#
# colvars molid top
# colvars configfile my_colvars.in
#
# start_colvar_display white;              # text color
# set_colvar_display position { 0 0 5 };   # change text position
# stop_colvar_display
#
###########################################################

proc start_colvar_display { {text_color blue} {cv all} {text_position {0 0 0}} {first none} {molid top} } {

  global text_pos text_col cv_handle

  set text_pos $text_position
  set text_col $text_color
  set cv_handle -1

  if {! [string compare $molid top]} {
    set molid [molinfo top]
  }

  # draw values right away
  update_colvar_display vmd_frame $molid w

  # set trace
  global vmd_frame 
  trace variable vmd_frame($molid) w update_colvar_display
  return
}


###########################################################3
proc stop_colvar_display { {molid top} } {
  global cv_handle

  if { $cv_handle > -1 } {
    graphics $molid delete $cv_handle
  }
  if {! [string compare $molid top]} {
    set molid [molinfo top]
  }
  trace vdelete vmd_frame($molid) w update_colvar_display
  return
}

###########################################################3
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


###########################################################3
proc update_colvar_display {name molid op} {
  # name == vmd_frame
  # molid == molecule id of the newly changed frame
  # op == w

  global text_pos text_col cv_handle
  set f [molinfo $molid get frame]

  colvars frame $f
  colvars update

  if { $cv_handle > -1 } {
    graphics $molid replace $cv_handle
  }

  set t ""
  foreach i [colvars printframe] {
    append t [format "%.2f " $i]
  }

  graphics $molid color $text_col
  set cv_handle [graphics $molid text $text_pos $t]
  return
}


