

proc ::cv_dashboard::createSettingsWindow { } {

  set w .cv_dashboard_window
  # set settings [toplevel $w.settings]
  # wm title $settings "Colvars Dashboard Settings"
  # wm protocol $settings WM_DELETE_WINDOW [list wm withdraw $settings]

  set settings $w.settings
  grid [frame $settings] -column 0 -columnspan 3 -sticky nsew

  set gridrow 0

  # Units
  incr gridrow
  grid [label $settings.unitTxt -text "Units:"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  ttk::combobox $settings.units -justify left -state readonly
  $settings.units configure -values [array names ::cv_dashboard::text_to_units]
  refresh_units
  grid $settings.units -row $gridrow -column 1 -columnspan 2 -pady 2 -padx 2 -sticky nsew
  bind $settings.units <<ComboboxSelected>> ::cv_dashboard::change_units

  incr gridrow
  grid [ttk::separator $settings.sep1 -orient horizontal] -row $gridrow -column 0 -columnspan 3 -pady 5 -sticky ew

  # Gradient display settings
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

  grid columnconfigure $settings 0 -weight 1
  grid columnconfigure $settings 1 -weight 1
  grid columnconfigure $settings 2 -weight 1

  grid remove $settings
  set ::cv_dashboard::settings_shown false
}


proc ::cv_dashboard::refresh_units {} {
  set settings .cv_dashboard_window.settings
  if [catch { set u [cv units] }] {
    # This catches cases where the module cannot be created because no molecule is loaded
    set u ""
  }
  set ::cv_dashboard::units $u
  if { $u == "" } {
    $settings.units set $::cv_dashboard::units_to_text(real)
  } else {
    $settings.units set $::cv_dashboard::units_to_text($u)
  }
}


# Change units if possible
proc ::cv_dashboard::change_units {} {
  set settings .cv_dashboard_window.settings
  set val [$settings.units get]
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
