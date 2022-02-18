

proc ::cv_dashboard::createSettingsWindow { } {

  set w .cv_dashboard_window

  set settings $w.tabs.settings
  grid [frame $settings] -column 0 -columnspan 3 -sticky nsew

  set gridrow 0

  # Graphics settings

  incr gridrow
  grid [ttk::separator $settings.sep1 -orient horizontal] -row $gridrow -column 0 -columnspan 3 -pady 5 -sticky ew

  incr gridrow
  grid [tk::label $settings.graphics_text -font $::cv_dashboard::font -text "Graphics settings"] -row $gridrow -column 0 -columnspan 3 -pady 2 -padx 2 -sticky nsew

  incr gridrow
  grid [tk::label $settings.material_text -font $::cv_dashboard::font -text "Display material:"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  ttk::combobox $settings.material -justify left -state readonly
  $settings.material configure -values [material list]
  grid $settings.material -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  $settings.material set "Opaque"

  incr gridrow
  grid [tk::label $settings.atom_radius_text -font $::cv_dashboard::font -text "Sphere radius:"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::entry $settings.atom_radius -textvariable ::cv_dashboard::atom_radius] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  set ::cv_dashboard::atom_radius 0.5

  incr gridrow
  grid [tk::label $settings.sel_text_label -font $::cv_dashboard::font -text "Intersect sel. text:"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::entry $settings.sel_text -textvariable ::cv_dashboard::sel_text] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  set ::cv_dashboard::sel_text "all"

  # Gradient display settings

  incr gridrow
  grid [ttk::separator $settings.sepgradients -orient horizontal] -row $gridrow -column 0 -columnspan 3 -pady 5 -sticky ew

  incr gridrow
  grid [tk::label $settings.gradients_text -font $::cv_dashboard::font -text "Gradients / Forces settings"] -row $gridrow -column 0 -columnspan 3 -pady 2 -padx 2 -sticky nsew

  incr gridrow
  grid [ttk::radiobutton $settings.grad_scale_choice_norm -text "Max. arrow length" -value "norm" -variable ::cv_dashboard::grad_scale_choice \
    -command {::cv_dashboard::update_shown_gradients; ::cv_dashboard::update_shown_forces} -padding "2 0 2 0"] \
    -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::entry $settings.grad_norm -textvariable ::cv_dashboard::grad_norm] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [tk::label $settings.grad_norm_unit -font $::cv_dashboard::font -text "Angstrom"] -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew
  bind $settings.grad_norm <<keyb_enter>> "$settings.grad_scale_choice_norm invoke; ::cv_dashboard::update_shown_gradients"

  incr gridrow
  grid [ttk::radiobutton $settings.grad_scale_choice_scale -text "Scaling factor" -value "scale" -variable ::cv_dashboard::grad_scale_choice \
    -command {::cv_dashboard::update_shown_gradients; ::cv_dashboard::update_shown_forces} -padding "2 0 2 0"] \
    -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::entry $settings.grad_scale -textvariable ::cv_dashboard::grad_scale] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [tk::label $settings.grad_scale_unit -font $::cv_dashboard::font -text "gradient: A*L/(cv/width)"] -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew
  bind $settings.grad_scale <<keyb_enter>> "$settings.grad_scale_choice_scale invoke; ::cv_dashboard::update_shown_gradients"

  $settings.grad_scale_choice_norm invoke ;# Default to norm

  incr gridrow
  grid [tk::label $settings.grad_radius_text -font $::cv_dashboard::font -text "Arrow radius:"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::entry $settings.grad_radius -textvariable ::cv_dashboard::grad_radius] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [tk::label $settings.grad_radius_unit -font $::cv_dashboard::font -text "Angstrom"] -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew
  set ::cv_dashboard::grad_radius 0.3

  # Rotation display settings

  incr gridrow
  grid [ttk::separator $settings.seprotations -orient horizontal] -row $gridrow -column 0 -columnspan 3 -pady 5 -sticky ew

  incr gridrow
  grid [tk::label $settings.rot_text -font $::cv_dashboard::font -text "Rotation display settings"] -row $gridrow -column 0 -columnspan 3 -pady 2 -padx 2 -sticky nsew

  incr gridrow
  grid [tk::label $settings.rot_scale_text -font $::cv_dashboard::font -text "Rotation object scale:"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::entry $settings.rot_scale -textvariable ::cv_dashboard::rot_scale] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  set ::cv_dashboard::rot_scale 1.0

  incr gridrow
  grid [tk::label $settings.rot_color_text -font $::cv_dashboard::font -text "Rotation object color:"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::combobox $settings.rot_color -justify left -state readonly] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  $settings.rot_color configure -values [colorinfo colors]
  $settings.rot_color set "yellow"

  # Volmap-specific settings

  incr gridrow
  grid [ttk::separator $settings.sepvolmap -orient horizontal] -row $gridrow -column 0 -columnspan 3 -pady 5 -sticky ew

  incr gridrow
  grid [tk::label $settings.volmap_text -font $::cv_dashboard::font -text "Volmaps settings"] -row $gridrow -column 0 -columnspan 3 -pady 2 -padx 2 -sticky nsew

  incr gridrow
  grid [tk::label $settings.volmap_contour_text -font $::cv_dashboard::font -text "Contour level:"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::entry $settings.volmap_contour -textvariable ::cv_dashboard::volmap_contour] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  set ::cv_dashboard::volmap_contour 0.5
  grid [tk::label $settings.volmap_contour_unit -font $::cv_dashboard::font -text "(% of min-max range)"] -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew

  incr gridrow
  grid [ttk::checkbutton $settings.volmap_periodic_x -text "+/-X images" -variable ::cv_dashboard::volmap_periodic_x] \
    -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::checkbutton $settings.volmap_periodic_y -text "+/-Y images" -variable ::cv_dashboard::volmap_periodic_y] \
    -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [ttk::checkbutton $settings.volmap_periodic_z -text "+/-Z images" -variable ::cv_dashboard::volmap_periodic_z] \
    -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew
  set ::cv_dashboard::volmap_periodic_x 0
  set ::cv_dashboard::volmap_periodic_y 0
  set ::cv_dashboard::volmap_periodic_z 0

  incr gridrow
  grid [ttk::separator $settings.sephist -orient horizontal] -row $gridrow -column 0 -columnspan 3 -pady 5 -sticky ew

  incr gridrow
  grid [tk::label $settings.hist_text -font $::cv_dashboard::font -text "Histogram settings"] -row $gridrow -column 0 -columnspan 3 -pady 2 -padx 2 -sticky nsew
  incr gridrow
  grid [tk::label $settings.nbins_text -font $::cv_dashboard::font -text "Number of bins:"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::entry $settings.nbins -textvariable ::cv_dashboard::nbins] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  set ::cv_dashboard::nbins 60

  grid columnconfigure $settings 0 -weight 1
  grid columnconfigure $settings 1 -weight 1
  grid columnconfigure $settings 2 -weight 1

  grid remove $settings
}


proc ::cv_dashboard::refresh_units {} {
  set main .cv_dashboard_window.tabs.main
  if [catch { set u [cv units] }] {
    # This catches cases where the module cannot be created because no molecule is loaded
    set u ""
  }
  set ::cv_dashboard::units $u
  if { $u == "" } {
    $main.units set $::cv_dashboard::units_to_text(real)
  } else {
    $main.units set $::cv_dashboard::units_to_text($u)
  }
}


# Change units if possible
proc ::cv_dashboard::change_units {} {
  set main .cv_dashboard_window.tabs.main
  set val [$main.units get]
  if {![info exists ::cv_dashboard::text_to_units($val)]} {
    puts "Bug error: trying to switch to unknown unit system $val"
    return
  }
  set new $::cv_dashboard::text_to_units($val)
  # Get up-to-date current setting
  refresh_units
  if {$new != $::cv_dashboard::units} {
    if {[run_cv list] != {}} {
      tk_messageBox -icon warning -title "Colvars Dashboard Warning"\
        -message "Warning: Changing units while colvars are defined.
Make sure the configuration of all variables is compatible with the new unit system. In particular, \
check any parameters with the dimension of a length or a force constant."
    }
  }
  cv units $new
  # Remember in global config
  dict set ::cv_dashboard::global_config units $new

  # Do a full reset to parse colvars again, in particular reading ref position files
  set cfg [get_whole_config]
  run_cv reset
  apply_config $cfg
}
