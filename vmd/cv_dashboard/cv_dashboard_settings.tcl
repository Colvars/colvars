

proc ::cv_dashboard::createSettingsWindow { } {

  set w .cv_dashboard_window

  set settings $w.tabs.settings
  grid [frame $settings] -column 0 -columnspan 3 -sticky nsew

  set gridrow 0

  # Graphics settings

  incr gridrow
  grid [ttk::separator $settings.sep1 -orient horizontal] -row $gridrow -column 0 -columnspan 3 -pady 5 -sticky ew

  incr gridrow
  grid [label $settings.graphics_text -text "Graphics settings"] -row $gridrow -column 0 -columnspan 3 -pady 2 -padx 2 -sticky nsew

  incr gridrow
  grid [label $settings.material_text -text "Display material:"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  ttk::combobox $settings.material -justify left -state readonly
  $settings.material configure -values [material list]
  grid $settings.material -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  $settings.material set "Opaque"

  incr gridrow
  grid [label $settings.atom_radius_text -text "Sphere radius:"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [tk::entry $settings.atom_radius -textvariable ::cv_dashboard::atom_radius] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  set ::cv_dashboard::atom_radius 0.5

  incr gridrow
  grid [label $settings.sel_text_label -text "Intersect sel. text:"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [tk::entry $settings.sel_text -textvariable ::cv_dashboard::sel_text] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  set ::cv_dashboard::sel_text "all"

  # Gradient display settings

  incr gridrow
  grid [ttk::separator $settings.sepgradients -orient horizontal] -row $gridrow -column 0 -columnspan 3 -pady 5 -sticky ew

  incr gridrow
  grid [label $settings.gradients_text -text "Gradients settings"] -row $gridrow -column 0 -columnspan 3 -pady 2 -padx 2 -sticky nsew

  incr gridrow
  grid [ttk::radiobutton $settings.grad_scale_choice_norm -text "Max. gradient norm" -value "norm" -variable ::cv_dashboard::grad_scale_choice \
    -command ::cv_dashboard::update_shown_gradients -padding "2 0 2 0"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [tk::entry $settings.grad_norm -textvariable ::cv_dashboard::grad_norm] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [label $settings.grad_norm_unit -text "Angstrom"] -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew
  bind $settings.grad_norm <<keyb_enter>> "$settings.grad_scale_choice_norm invoke; ::cv_dashboard::update_shown_gradients"

  incr gridrow
  grid [ttk::radiobutton $settings.grad_scale_choice_scale -text "Grad. scaling factor" -value "scale" -variable ::cv_dashboard::grad_scale_choice \
    -command ::cv_dashboard::update_shown_gradients -padding "2 0 2 0"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [tk::entry $settings.grad_scale -textvariable ::cv_dashboard::grad_scale] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [label $settings.grad_scale_unit -text "A * L / (cv / width)"] -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew
  bind $settings.grad_scale <<keyb_enter>> "$settings.grad_scale_choice_scale invoke; ::cv_dashboard::update_shown_gradients"

  $settings.grad_scale_choice_norm invoke ;# Default to norm

  incr gridrow
  grid [label $settings.grad_radius_text -text "Arrow radius:"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [tk::entry $settings.grad_radius -textvariable ::cv_dashboard::grad_radius] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [label $settings.grad_radius_unit -text "Angstrom"] -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew
  set ::cv_dashboard::grad_radius 0.3

  # Rotation display settings

  incr gridrow
  grid [ttk::separator $settings.seprotations -orient horizontal] -row $gridrow -column 0 -columnspan 3 -pady 5 -sticky ew

  incr gridrow
  grid [label $settings.rot_text -text "Rotation display settings"] -row $gridrow -column 0 -columnspan 3 -pady 2 -padx 2 -sticky nsew

  incr gridrow
  grid [label $settings.rot_scale_text -text "Rotation object scale:"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [tk::entry $settings.rot_scale -textvariable ::cv_dashboard::rot_scale] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  set ::cv_dashboard::rot_scale 1.0

  incr gridrow
  grid [label $settings.rot_color_text -text "Rotation object color:"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::combobox $settings.rot_color -justify left -state readonly] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  $settings.rot_color configure -values [colorinfo colors]
  $settings.rot_color set "yellow"

  # Volmap-specific settings

  incr gridrow
  grid [ttk::separator $settings.sepvolmap -orient horizontal] -row $gridrow -column 0 -columnspan 3 -pady 5 -sticky ew

  incr gridrow
  grid [label $settings.volmap_text -text "Volmaps settings"] -row $gridrow -column 0 -columnspan 3 -pady 2 -padx 2 -sticky nsew

  incr gridrow
  grid [label $settings.volmap_contour_text -text "Contour level:"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [tk::entry $settings.volmap_contour -textvariable ::cv_dashboard::volmap_contour] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  set ::cv_dashboard::volmap_contour 0.5
  grid [label $settings.volmap_contour_unit -text "(% of min-max range)"] -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew

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

  grid columnconfigure $settings 0 -weight 1
  grid columnconfigure $settings 1 -weight 1
  grid columnconfigure $settings 2 -weight 1

  grid remove $settings
  set ::cv_dashboard::settings_shown false
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
    if {[run_cv list] == {}} {
      cv units $new
    } else {
      tk_messageBox -icon error -title "Colvars Dashboard Error"\
        -message "Cannot change units while colvars are defined.
You can either:
1) delete all colvars, or
2) edit their configuration (<Ctrl-a> <e>), checking all parameters \
for consistency with desired units, and add the following line:
units $new"
      return
    }
  }
  # Refresh Combo box
  refresh_units
  refresh_values
}
