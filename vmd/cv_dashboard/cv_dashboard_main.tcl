#################################################################
# Main UI: the dashboard
#################################################################

# -*- tcl-indent-level: 2; -*-


# Create main window
proc ::cv_dashboard::createWindow {} {

  set w [toplevel .cv_dashboard_window]
  wm title $w "Colvars dashboard"
  wm protocol $w WM_DELETE_WINDOW { ::cv_dashboard::quit }

  # Top bars of buttons
  set gridrow 0
  grid [ttk::button $w.helpB -text "Online Help" -command {::cv_dashboard::invokeBrowser "http://colvars.github.io/colvars-refman-vmd/colvars-refman-vmd.html#sec:dashboard"} -padding "2 0 2 0"] \
    -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $w.aboutB -text "About" -command ::cv_dashboard::about -padding "2 0 2 0"] \
    -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $w.quit -text "Quit" -command ::cv_dashboard::quit -padding "2 0 2 0"] \
    -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew

  incr gridrow
  grid [ttk::button $w.load -text "Load" -command ::cv_dashboard::load -padding "2 0 2 0"] \
    -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $w.save -text "Save" -command ::cv_dashboard::save -padding "2 0 2 0"] \
    -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $w.reset -text "Reset" -command ::cv_dashboard::reset -padding "2 0 2 0"] \
    -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew

  # Table of colvars
  ttk::treeview $w.cvtable -selectmode extended -show {headings tree} -height 8
  $w.cvtable configure -column val
  $w.cvtable column #0 -width 50 -stretch 1 -anchor w
  $w.cvtable column val -width 150 -stretch 1 -anchor w

  $w.cvtable heading #0 -text "colvar name"
  $w.cvtable heading val -text "value"

  bind $w.cvtable <Button-3> {::cv_dashboard::cvContextMenu %x %y %X %Y}
  bind $w.cvtable <Button-1> {::cv_dashboard::tableClicked cvtable %x %y}

  bind $w.cvtable <Control-e> ::cv_dashboard::edit_cv
  bind $w.cvtable <Control-a> { .cv_dashboard_window.cvtable selection set $::cv_dashboard::cvs }
  bind $w.cvtable <Control-n> ::cv_dashboard::add_cv

  event add <<keyb_enter>> <Return>   ;# Combine Return and keypad-Enter into a single virtual event
  event add <<keyb_enter>> <KP_Enter>

  if { [info patchlevel] != "8.5.6" } {
    $w.cvtable tag configure parity0 -background white
    $w.cvtable tag configure parity1 -background grey94
  }

  incr gridrow
  grid $w.cvtable -row $gridrow -column 0 -sticky news -columnspan 3 -pady 10
  # The colvar table expands and shrinks with the window height
  grid rowconfigure $w $gridrow -weight 1 -minsize 20

  # Tabs
  incr gridrow
  grid [ttk::notebook $w.tabs] -row $gridrow -column 0 -columnspan 3 -sticky news -padx 0
  grid rowconfigure $w $gridrow -weight 1 -minsize 20
  tk::frame $w.tabs.main

  # Editing
  set gridrow 0  ;# New grid
  set main .cv_dashboard_window.tabs.main

  grid [ttk::separator $main.sep1 -orient horizontal] -row $gridrow -column 0 -columnspan 3 -pady 5 -sticky ew

  incr gridrow
  grid [label $main.actions_text -text "Colvar list actions"] -row $gridrow -column 0 -columnspan 3 -pady 2 -padx 2 -sticky nsew

  incr gridrow

  grid [ttk::button $main.edit -text "Edit \[Ctrl-e\]" -command ::cv_dashboard::edit_cv -padding "2 0 2 0"] \
    -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $main.add -text "New \[Ctrl-n\]" -command ::cv_dashboard::add_cv -padding "2 0 2 0"] \
    -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $main.del -text "Delete" -command ::cv_dashboard::del_cv -padding "2 0 2 0"] \
    -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew
  incr gridrow
  grid [ttk::button $main.refresh -text "Refresh list \[F5\]" -command ::cv_dashboard::refresh_table -padding "2 0 2 0"] -row $gridrow -column 0 -columnspan 3 -pady 2 -padx 2 -sticky nsew

  # Plots
  incr gridrow
  grid [ttk::separator $main.sep_plots -orient horizontal] -row $gridrow -column 0 -columnspan 3 -pady 5 -sticky ew
  incr gridrow
  grid [label $main.viz_text -text "Plots and real-time visualizations"] -row $gridrow -column 0 -columnspan 3 -pady 2 -padx 2 -sticky nsew
  incr gridrow
  grid [ttk::button $main.plot -text "Timeline plot" -command ::cv_dashboard::plot -padding "2 0 2 0"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $main.plot2cv -text "Pairwise plot" -command {::cv_dashboard::plot 2cv} -padding "2 0 2 0"] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew

  user add key F5 ::cv_dashboard::refresh_table
  bind $main <F5> ::cv_dashboard::refresh_table

  # Add trajectory animation bindings to Dashboard and VMD's OpenGL window
  traj_animation_bindings $w
  user add key Left  { ::cv_dashboard::chg_frame -1 }
  user add key Right { ::cv_dashboard::chg_frame 1 }
  user add key Home  { ::cv_dashboard::chg_frame start }
  user add key End   { ::cv_dashboard::chg_frame end }

  # Atom group display
  incr gridrow
  grid [ttk::button $main.show_atoms -text "Show atoms" -command {::cv_dashboard::show_atoms_selected} -padding "2 0 2 0"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $main.hide_atoms -text "Hide atoms" -command {::cv_dashboard::hide_atoms_selected} -padding "2 0 2 0"] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $main.hide_all_atoms -text "Hide all atoms" -command {::cv_dashboard::hide_all_atoms} -padding "2 0 2 0"] -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew

  # Gradient display
  incr gridrow
  grid [ttk::button $main.show_gradients -text "Show gradients" -command {::cv_dashboard::show_gradients [::cv_dashboard::selected_colvars]} \
    -padding "2 0 2 0"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $main.hide_gradients -text "Hide gradients" -command {::cv_dashboard::hide_gradients} -padding "2 0 2 0"] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $main.hide_all_gradients -text "Hide all grads" -command {::cv_dashboard::hide_all_gradients} -padding "2 0 2 0"] -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew

  # Create and hide rotation menu (shows itself as needed when refreshing the colvar table)
  incr gridrow
  createRotationMenu $gridrow

  # Create and hide volmap menu (shows itself as needed when refreshing the colvar table)
  incr gridrow
  createVolmapMenu $gridrow

  # Wizards
  incr gridrow
  grid [ttk::separator $main.sep_auto -orient horizontal] -row $gridrow -column 0 -columnspan 3 -pady 5 -sticky ew
  incr gridrow
  grid [label $main.auto_text -text "Automatic colvars"] -row $gridrow -column 0 -columnspan 3 -pady 2 -padx 2 -sticky nsew
  incr gridrow
  grid [ttk::button $main.cvfromprotein -text "Protein auto-colvars" -command ::cv_dashboard::protein_cvs -padding "2 0 2 0"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $main.cvfromlabels -text "Colvars from VMD labels" -command ::cv_dashboard::cvs_from_labels -padding "2 0 2 0"] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew

  # General options
  incr gridrow
  grid [ttk::separator $main.sep_options -orient horizontal] -row $gridrow -column 0 -columnspan 3 -pady 5 -sticky ew
  incr gridrow
  grid [label $main.options_text -text "General options"] -row $gridrow -column 0 -columnspan 3 -pady 2 -padx 2 -sticky nsew

  incr gridrow
  grid [label $main.molTxt -text "Molecule:"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  ttk::combobox $main.mol -justify left -state readonly
  $main.mol configure -values [molinfo list]
  if { $::cv_dashboard::mol != -1 } {
    $main.mol set $::cv_dashboard::mol
  }
  trace add variable ::vmd_initialize_structure write ::cv_dashboard::update_mol_list
  grid $main.mol -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  bind $main.mol <<ComboboxSelected>> ::cv_dashboard::change_mol
  grid [ttk::button $main.molTop -text "Switch to top" -command {::cv_dashboard::switch_to_top_mol} -padding "2 0 2 0"] -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew


  # Units
  incr gridrow
  grid [label $main.unitTxt -text "Units:"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  ttk::combobox $main.units -justify left -state readonly
  $main.units configure -values [array names ::cv_dashboard::text_to_units]
  refresh_units
  grid $main.units -row $gridrow -column 1 -columnspan 2 -pady 2 -padx 2 -sticky nsew
  bind $main.units <<ComboboxSelected>> ::cv_dashboard::change_units

  # Frame display and track checkbox
  incr gridrow
  grid [label $main.frameTxt -text "Frame:"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [label $main.frame -textvariable ::cv_dashboard::current_frame] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [ttk::checkbutton $main.trackFrame -text "Track VMD frame" -command ::cv_dashboard::change_track_frame -variable ::cv_dashboard::track_frame] \
    -row $gridrow -column 2  -pady 2 -padx 2 -sticky nsew

  # Create and hide Settings window to create all associated variables
  createSettingsWindow

  # Create and hide Biases tab
  createBiasesTab

  $w.tabs add $w.tabs.main -text "Actions" -sticky news
  $w.tabs add $w.tabs.biases -text "Biases" -sticky news
  $w.tabs add $w.tabs.settings -text "Settings" -sticky news
  # $w.tabs add $w.tabs.stats -text "Force/energy stats" -sticky news

  grid columnconfigure $main 0 -weight 1
  grid columnconfigure $main 1 -weight 1
  grid columnconfigure $main 2 -weight 1

  grid columnconfigure $w 0 -weight 1
  grid columnconfigure $w 1 -weight 1
  grid columnconfigure $w 2 -weight 1

  refresh_table
  change_track_frame ;# activate tracking if necessary

  return $w
}

proc ::cv_dashboard::createBiasesTab {} {

  set biases .cv_dashboard_window.tabs.biases

  set gridrow 0
  grid [frame $biases] -row $gridrow -column 0 -columnspan 3 -sticky nsew

  # Table of biases
  ttk::treeview $biases.bias_table -selectmode extended -show {headings tree} -height 3
  $biases.bias_table configure -columns { val colvars }
  $biases.bias_table column #0 -width 40 -stretch 1 -anchor w
  $biases.bias_table column val -width 20 -stretch 1 -anchor w
  $biases.bias_table column colvars -width 60 -stretch 1 -anchor w

  $biases.bias_table heading #0 -text "bias name"
  $biases.bias_table heading val -text "energy"
  $biases.bias_table heading colvars -text "colvars"

  bind $biases.bias_table <Control-a> { .cv_dashboard_window.tabs.biases.bias_table selection set $::cv_dashboard::biases }
  bind $biases.bias_table <Button-1> {::cv_dashboard::tableClicked tabs.biases.bias_table %x %y}
  bind $biases.bias_table <Button-3> {::cv_dashboard::biasContextMenu %x %y %X %Y}

  event add <<keyb_enter>> <Return>   ;# Combine Return and keypad-Enter into a single virtual event
  event add <<keyb_enter>> <KP_Enter>

  if { [info patchlevel] != "8.5.6" } {
    $biases.bias_table tag configure parity0 -background white
    $biases.bias_table tag configure parity1 -background grey94
  }

  incr gridrow
  grid $biases.bias_table -row $gridrow -column 0 -sticky news -columnspan 3
  grid rowconfigure $biases $gridrow -weight 1 -minsize 20

  incr gridrow
  grid [ttk::separator $biases.sep1 -orient horizontal] -row $gridrow -column 0 -columnspan 3 -pady 5 -sticky ew

  incr gridrow
  grid [label $biases.actions_text -text "Bias list actions"] -row $gridrow -column 0 -columnspan 3 -pady 2 -padx 2 -sticky nsew

  incr gridrow
  grid [ttk::button $biases.edit -text "Edit bias" -command ::cv_dashboard::edit_bias -padding "2 0 2 0"] \
    -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $biases.add -text "New bias" -command ::cv_dashboard::add_bias -padding "2 0 2 0"] \
    -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $biases.del -text "Delete bias" -command ::cv_dashboard::del_bias -padding "2 0 2 0"] \
    -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew

  incr gridrow
  grid [ttk::button $biases.plot -text "Energy timeline" -command ::cv_dashboard::plot_bias_energy -padding "2 0 2 0"] \
    -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $biases.refresh -text "Refresh list" -command ::cv_dashboard::refresh_bias_table -padding "2 0 2 0"] \
    -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew

  if { [string compare [run_cv version] "2021-12-20"] >= 0 } {
    # Force display
    incr gridrow
    grid [ttk::button $biases.show_forces -text "Show bias forces" -command {::cv_dashboard::show_forces [::cv_dashboard::selected_biases]} \
      -padding "2 0 2 0"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
    grid [ttk::button $biases.hide_forces -text "Hide bias forces" -command {::cv_dashboard::hide_forces} -padding "2 0 2 0"] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
    grid [ttk::button $biases.hide_all_forces -text "Hide all forces" -command {::cv_dashboard::hide_all_forces} -padding "2 0 2 0"] -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew
  }
  # Stats "tab" is now included in biases tab, appending to current grid
  createStatsTab $gridrow

  grid columnconfigure $biases 0 -weight 1
  grid columnconfigure $biases 1 -weight 1
  grid columnconfigure $biases 2 -weight 1

  grid remove $biases
}


# Display context menu about specific bias
# Takes coordinates within widget and within window
proc ::cv_dashboard::biasContextMenu { x y wX wY } {
  set w .cv_dashboard_window.tabs.biases
  set menu $w.biasMenu

  set biases [selected_biases]

  # Add any bias under mouse cursor
  set b [$w.bias_table identify row $x $y]
  if { [llength $b] == 1 && [lsearch $biases $b] == -1 } {
    lappend biases $b
  }

  if { [winfo exists $menu] } {
    destroy $menu
  }
  menu $menu -tearoff 0

  if { [llength $biases] == 0 } {
    $menu add command -label New -command ::cv_dashboard::add_bias
  } else {
    $menu add command -label Edit -command [list ::cv_dashboard::edit_bias false $biases]
    $menu add command -label Delete -command [list ::cv_dashboard::del_bias $biases]
  }
  tk_popup $menu $wX $wY
}


proc ::cv_dashboard::createStatsTab { gridrow } {

  # Merge stats tab with biases tab
  set stats .cv_dashboard_window.tabs.biases

  incr gridrow
  grid [ttk::separator $stats.sep_stats -orient horizontal] -row $gridrow -column 0 -columnspan 3 -pady 5 -sticky ew

  incr gridrow
  grid [label $stats.stats_title -text "Energy and force statistics"] -row $gridrow -column 0 -columnspan 3 -pady 2 -padx 2 -sticky nsew

  incr gridrow
  # Energy/Force display
  set ::cv_dashboard::colvar_energy 0.0
  grid [label $stats.energyTxt -text "Total energy:"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [tk::entry $stats.energy -textvariable ::cv_dashboard::colvar_energy -state readonly] \
    -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew

  incr gridrow
  set ::cv_dashboard::atom_forces_rms 0.0
  set ::cv_dashboard::atom_forces_max 0.0
  set ::cv_dashboard::atom_forces_max_id -1
  grid [label $stats.rmsForceTxt -text "RMS force:"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [tk::entry $stats.rmsForce -textvariable ::cv_dashboard::atom_forces_rms -state readonly] \
    -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [label $stats.maxForceIDTxt -text "Max force atom index:"] -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew

  incr gridrow
  grid [label $stats.maxForceTxt -text "Max force:"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [tk::entry $stats.maxForce -textvariable ::cv_dashboard::atom_forces_max -state readonly] \
    -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [tk::entry $stats.maxForceID -textvariable ::cv_dashboard::atom_forces_max_id -state readonly] \
    -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew
}

# Open or close the volmap sub-panel
proc ::cv_dashboard::toggleVolmapMenu {} {
  set w .cv_dashboard_window.tabs.main

  set volmaps false
  foreach cv $::cv_dashboard::cvs {
    if { [is_volmap $cv] } {
      set volmaps true
      break
    }
  }

  if { $volmaps } {
    grid $w.show_volmaps $w.hide_volmaps $w.hide_all_volmaps
  } else {
    grid remove $w.show_volmaps $w.hide_volmaps $w.hide_all_volmaps
  }
}


# Open or close the volmap sub-panel
proc ::cv_dashboard::toggleRotationMenu {} {
  set w .cv_dashboard_window.tabs.main

  set rotations false
  foreach cv $::cv_dashboard::cvs {
    if { [is_unit_quaternion $cv] } {
      set rotations true
      break
    }
  }

  if { $rotations } {
    grid $w.show_rotation $w.hide_rotation
  } else {
    grid remove $w.show_rotation $w.hide_rotation
  }
}


# Display context menu about specific colvar
# Takes coordinates within widget and within window
proc ::cv_dashboard::cvContextMenu { x y wX wY } {
  set w .cv_dashboard_window
  set menu $w.cvMenu

  set cvs [selected_colvars]

  # Add any colvar under mouse cursor
  # reduce scalar components to their parent vector CV
  set cv [lindex [$w.cvtable identify row $x $y] 0]
  if { [llength $cv] == 1 && [lsearch $cvs $cv] == -1 } {
    lappend cvs $cv
  }

  set volmaps [list]
  set rotations [list]
  foreach cv $cvs {
    if { [is_volmap $cv] } {
      lappend volmaps $cv
    }
    if { [is_unit_quaternion $cv] } {
      lappend rotations $cv
    }
  }

  if { [winfo exists $menu] } {
    destroy $menu
  }
  menu $menu -tearoff 0

  if { [llength $cvs] == 0 } {
    $menu add command -label New -command ::cv_dashboard::add_cv
  } else {
    if { [llength $volmaps] > 0 } {
      $menu add command -label "Show volmap" -command [list ::cv_dashboard::show_volmaps $volmaps]
      $menu add command -label "Hide volmap" -command [list ::cv_dashboard::hide_volmaps $volmaps]
    }
    if { [llength $rotations] > 0 } {
      $menu add command -label "Show rotation" -command [list ::cv_dashboard::start_rotation_display $rotations]
      $menu add command -label "Hide rotation" -command [list ::cv_dashboard::stop_rotation_display]
    }
    $menu add command -label Edit -command [list ::cv_dashboard::edit_cv false $cvs]
    $menu add command -label Delete -command [list ::cv_dashboard::del_cv $cvs]
    $menu add command -label "Add harmonic bias" -command [list ::cv_dashboard::add_bias $cvs]
  }
  tk_popup $menu $wX $wY
}


# Takes coordinates within widget
proc ::cv_dashboard::tableClicked { table x y } {
  set t .cv_dashboard_window.$table

  # colvar / bias under mouse
  set item [$t identify row $x $y]

  # Deselect all if clicked on nothing
  if { [llength $item] == 0 } {
    $t selection set [list]
  }
}


# Display help window with basic information
proc ::cv_dashboard::about {} {
  help_window .cv_dashboard_window "About the Colvars Dashboard" "About the Colvars Dashboard" \
"Colvars Dashboard, version [package require cv_dashboard]

Based on the Colvars Module version [run_cv version]
In [vmdinfo versionmsg]
Running Tcl/Tk [info patchlevel]

Jérôme Hénin (henin@ibpc.fr) and the Colvars developers.

J. Hénin, L. J. S. Lopes, and G. Fiorin. Human learning for molecular simulations: the Collective Variables Dashboard in VMD.
http://arxiv.org/abs/2110.08758
"
}


proc ::cv_dashboard::quit {} {
  # remove the trace we put in place, so they don't accumulate
  # if loaded multiple times
  catch {
    trace remove variable ::vmd_frame($::cv_dashboard::mol) write ::cv_dashboard::update_frame
  }
  wm withdraw .cv_dashboard_window
}


# Refresh the table with a list of existing CVs and their values
proc ::cv_dashboard::refresh_table {} {

  set w .cv_dashboard_window

  foreach i [$w.cvtable children {}] {
    $w.cvtable delete $i
  }

  if [catch { set ::cv_dashboard::cvs [cv list]}] {
    # We were unable to fetch the list of colvars
    # CVM is probably not initialized or there is no molecule loaded
    set ::cv_dashboard::cvs {}
    return
  }

  # Get fresh coordinates from VMD
  run_cv update

  set parity 1
  foreach cv $::cv_dashboard::cvs {
    # tag odd and even rows for alternating background colors
    set parity [expr 1-$parity]
    $w.cvtable insert {} end -id $cv -text $cv -tag parity$parity

    set val [run_cv colvar $cv update]
    $w.cvtable set $cv val [format_value $val]

    # Add sub-elements for vector colvars
    set size [llength $val]
    if { $size > 1 } {
      for {set i 1} {$i <= $size} {incr i} {
        set n "$cv $i"
        $w.cvtable insert $cv end -id $n -text $n -tag parity$parity
        $w.cvtable set $n val [format_value [lindex $val [expr $i-1]]]
      }
    }
  }

  refresh_energy_forces
  refresh_bias_table

  toggleRotationMenu
  toggleVolmapMenu

  update_frame internal $::cv_dashboard::mol w
}


# Refresh the table with a list of colvar values
proc ::cv_dashboard::refresh_values {} {
  run_cv update
  set w .cv_dashboard_window

  foreach cv [$w.cvtable children {}] {
    set val [run_cv colvar $cv update]
    $w.cvtable set $cv val [format_value $val]

    set size [llength $val]
    if { $size > 1 } {
      for {set i 1} {$i <= $size} {incr i} {
        set n "$cv $i"
        $w.cvtable set $n val [format_value [lindex $val [expr $i-1]]]
      }
    }
  }

  # TODO refresh only if tab is open and upon opening the tab
  set biases .cv_dashboard_window.tabs.biases
  foreach bias [$biases.bias_table children {}] {
    set val [run_cv bias $bias update]
    $biases.bias_table set $bias val [format_value $val]
  }

  refresh_energy_forces
}


# Refresh the table with a list of existing biases
proc ::cv_dashboard::refresh_bias_table {} {

  set biases .cv_dashboard_window.tabs.biases

  foreach i [$biases.bias_table children {}] {
    $biases.bias_table delete $i
  }
  set ::cv_dashboard::biases [cv list biases]

  set parity 1
  foreach bias $::cv_dashboard::biases {
    # tag odd and even rows for alternating background colors
    set parity [expr 1-$parity]
    $biases.bias_table insert {} end -id $bias -text $bias -tag parity$parity

    set val [run_cv bias $bias update]
    $biases.bias_table set $bias val [format_value $val]
    set cfg [run_cv bias $bias getconfig]
    set cvs ""
    regexp -line -nocase {^\s*colvars\s+(.*)} $cfg match cvs
    $biases.bias_table set $bias colvars $cvs
  }
}


proc ::cv_dashboard::refresh_energy_forces {} {
  catch {set ::cv_dashboard::colvar_energy [cv getenergy]}

  if { [string compare [run_cv version] "2021-03-02"] >= 0 } {
    set ::cv_dashboard::atom_forces_rms [run_cv getatomappliedforcesrms]
    set ::cv_dashboard::atom_forces_max [run_cv getatomappliedforcesmax]
    set ::cv_dashboard::atom_forces_max_id [run_cv getatomappliedforcesmaxid]
  } else {
    set ::cv_dashboard::atom_forces_rms "n/a"
    set ::cv_dashboard::atom_forces_max "n/a"
    set ::cv_dashboard::atom_forces_max_id "n/a"
  }
}


# Format colvar values returned by cv for display in table
proc ::cv_dashboard::format_value val {
  if {[llength $val] == 1} {
    if [catch { set str [format "%.4g" $val] }] {
      # values could not be formatted, maybe it is "nan"? Show as is.
      return $val
    } else {
      return $str
    }
  } else {
    set s "("
    foreach v [lrange $val 0 end-1] {
      append s "[format_value $v], "
    }
    append s "[format_value [lindex $val end]])"
    return $s
  }
}


# Return list of selected colvars in the table, or the lone colvar if there is just one
proc ::cv_dashboard::selected_colvars {} {
  set w .cv_dashboard_window

  if { [llength $::cv_dashboard::cvs] == 1 } {
    return [lindex $::cv_dashboard::cvs 0]
  } else {
    set cvs {}
    set prev ""
    foreach i [$w.cvtable selection] {
      # reduce scalar components to their parent vector CV
      set cv [lindex $i 0]
      # Skip series of duplicates
      if {$cv != $prev} {
        lappend cvs $cv
        set prev $cv
      }
    }
    return $cvs
  }
}


# Return selected scalar components for a given selected colvar
proc ::cv_dashboard::selected_comps { cv } {
  set w .cv_dashboard_window
  set comps {}

  foreach i [$w.cvtable selection] {
    # get comp index of matching cv elements
    if {[lindex $i 0] == $cv && [llength $i] == 2 } {
      lappend comps [expr [lindex $i 1] - 1]
    }
  }

  # If none selected, select all
  if { [llength $comps] == 0 } {
    set val [run_cv colvar $cv value]
    set i 0
    foreach v $val { lappend comps $i; incr i }
  }
  # Order of selected components is not guaranteed
  return [lsort $comps]
}


# Enable or disable real-time tracking of VMD frame
proc ::cv_dashboard::change_track_frame {} {
  set molid $::cv_dashboard::mol
  if { $molid == -1 } { return }
  if {$::cv_dashboard::track_frame} {
    trace add variable ::vmd_frame($molid) write ::cv_dashboard::update_frame
    update_frame internal $molid w
  } else {
    trace remove variable ::vmd_frame($molid) write ::cv_dashboard::update_frame
  }
}


# Load config from file
proc ::cv_dashboard::load {} {
  if { [info exists ::cv_dashboard::config_dir] } {
    set path [tk_getOpenFile -filetypes {{"Colvars cfg" .colvars} {"Colvars cfg" .in} {"Gromacs Colvars cfg" .dat} {"All files" *}} \
        -initialdir $::cv_dashboard::config_dir]
  } else {
    set path [tk_getOpenFile -filetypes {{"Colvars cfg" .colvars} {"Colvars cfg" .in} {"Gromacs Colvars cfg" .dat} {"All files" *}} \
        -initialdir [pwd]]
  }
  if [string compare $path ""] {
    set in [open $path r]
    set cfg [read -nonewline $in]
    close $in
    apply_config $cfg
    # Save directory for next invocation of this dialog
    set ::cv_dashboard::config_dir [file dirname $path]
  }
}


# Save config string of whole Colvars module to file
proc ::cv_dashboard::save {} {
  if { [info exists ::cv_dashboard::config_dir] } {
    set initialdir $::cv_dashboard::config_dir
  } else {
    set initialdir [pwd]
  }
  set path [tk_getSaveFile -filetypes {{"Colvars cfg" .colvars} {"Colvars cfg" .in} {"Gromacs Colvars cfg" .dat} {"All files" *}} \
        -initialdir $initialdir]

  if [string compare $path ""] { ;# Empty if the dialog is closed without confirmation
    # Save directory for next invocation of this dialog
    set ::cv_dashboard::config_dir [file dirname $path]

    set o [open $path w]
    puts $o [get_whole_config]
    close $o
  }
}


# Delete currently selected colvars
proc ::cv_dashboard::del_cv { {cvs "" } } {
  if { [llength $cvs] < 1 } {
    set cvs [selected_colvars]
  }
  foreach cv $cvs {
    # workaround bug in colvars pre-2018-07-02
    if {[string compare [run_cv version] "2018-07-02"] < 0} {
      foreach b [run_cv list biases] {
        if [string match *colvars*$c* [run_cv bias $b getconfig]] {
          puts "Bug workaround: deleting bias $b"
          run_cv bias $b delete
        }
      }
    }
    run_cv colvar $cv delete
  }
  refresh_table
}


# Reset cvm: hard delete
proc ::cv_dashboard::reset {} {
  set molid $::cv_dashboard::mol
  if { $::vmd_initialize_structure($molid) == 1 } {
    # If our molecule still exists,
    # remove all graphical objects which would be orphaned
    ::cv_dashboard::hide_all_atoms
    ::cv_dashboard::hide_all_gradients
    ::cv_dashboard::hide_all_forces
  }

  # Run cv delete silently to be less verbose if module was already deleted
  catch { cv delete }
  if { $molid != -1 } {
    run_cv molid $molid
  }
  set ::cv_dashboard::colvar_configs [dict create]
  set ::cv_dashboard::bias_configs [dict create]
  set ::cv_dashboard::global_config [dict create]
  set ::cv_dashboard::global_comments ""
  refresh_table
  refresh_units
}


# Return list of selected biases in the table
proc ::cv_dashboard::selected_biases {} {

  set biases .cv_dashboard_window.tabs.biases

  return [$biases.bias_table selection]
}


# Delete currently selected biases
proc ::cv_dashboard::del_bias { {biases "" } } {
  if { [llength $biases] < 1 } {
    set biases [selected_biases]
  }

  foreach bias $biases {
    run_cv bias $bias delete
  }
  refresh_bias_table
}


#################################################################
# Display atoms
#################################################################


# Display atoms in groups for selected colvars
proc ::cv_dashboard::show_atoms_selected {} {
  show_atoms [selected_colvars]
}


# Display atoms in groups for selected colvars
proc ::cv_dashboard::show_atoms { colvars } {
  set w .cv_dashboard_window
  set color 0
  set ci 0
  foreach cv $colvars {
    if { [info exists ::cv_dashboard::atom_rep($cv)]} {
      # Refreshing for this colvar: let's delete and re-create
      hide_atoms $cv
    }
    incr ci
    set all_groups [run_cv colvar $cv getatomgroups]
    # Remove characters such as <> which are parsed as special in VMD selection texts
    set sanitized_cvname [regsub -all {[^a-zA-Z0-9_@]} $cv {} ]
    set i 0
    set macros {}
    set repnames {}
    foreach list $all_groups {
      incr i
      # dummyAtoms will return empty lists
      if {[llength $list] > 0} {
        set group "${sanitized_cvname}_group_${i}"
        # resolve ambiguous names due to colvar name sanitization
        while {[lsearch [array names ::cv_dashboard::atom_rep] $group] > -1} {
          append sanitized_cvname "_"
          set group "${sanitized_cvname}_group_${i}"
        }
        lappend macros $group
        atomselect macro $group "($::cv_dashboard::sel_text) and (index $list)"
        mol color ColorID $color
        mol representation VDW [$w.tabs.settings.atom_radius get] 12
        mol selection "$group"
        mol material [$w.tabs.settings.material get]
        mol addrep $::cv_dashboard::mol
        set repid [expr [molinfo $::cv_dashboard::mol get numreps] - 1]
        lappend repnames [mol repname $::cv_dashboard::mol $repid]
        incr color
      }
    }
    set ::cv_dashboard::atom_rep($cv) [list $macros $repnames]
  }
}


# Hide atoms for selected colvars
proc ::cv_dashboard::hide_atoms_selected {} {
  hide_atoms [selected_colvars]
}


# Hide atoms for colvars given as parameters
proc ::cv_dashboard::hide_atoms { colvars } {
  foreach cv $colvars {
    if { [info exists ::cv_dashboard::atom_rep($cv)] } {
      lassign $::cv_dashboard::atom_rep($cv) macros repnames
      foreach m $macros {
        atomselect delmacro $m
      }
      foreach r $repnames {
        mol delrep [mol repindex $::cv_dashboard::mol $r] $::cv_dashboard::mol
      }
      unset ::cv_dashboard::atom_rep($cv)
    }
  }
}


# Remove all atom representations
proc ::cv_dashboard::hide_all_atoms {} {
  foreach {cv data} [array get ::cv_dashboard::atom_rep] {
    lassign $data macros repnames
    foreach m $macros {
      atomselect delmacro $m
    }
    foreach r $repnames {
      mol delrep [mol repindex $::cv_dashboard::mol $r] $::cv_dashboard::mol
    }
  }
  array unset ::cv_dashboard::atom_rep *
}


#################################################################
# Display volmaps
#################################################################


# Display atoms in groups for selected colvars
proc ::cv_dashboard::show_volmaps_selected {} {
  show_volmaps [selected_colvars]
}


# Display volmaps used by selected colvars
proc ::cv_dashboard::show_volmaps { colvars } {
  set w .cv_dashboard_window
  set color 0
  foreach cv $colvars {
    if { [info exists ::cv_dashboard::volmap_rep($cv)]} {
      hide_volmaps $cv
    }
    catch {set cv_volmaps [cv colvar $cv getvolmapids]}
    set repnames {}
    foreach volid $cv_volmaps {
      if { $volid >= 0 } {
        set threshold [expr $::cv_dashboard::volmap_contour * [vecsum \
          [voltool info minmax -mol $::cv_dashboard::mol -vol ${volid}]]]
        mol color ColorID $color
        mol selection all  ;# Must provide some selection text
        mol representation Isosurface ${threshold} ${volid} 2 0 0 1
        mol material [$w.tabs.settings.material get]
        mol addrep $::cv_dashboard::mol
        set repid [expr [molinfo $::cv_dashboard::mol get numreps] - 1]
        set periodic_string ""
        if { $::cv_dashboard::volmap_periodic_x } {
          set periodic_string ${periodic_string}"xX"
        }
        if { $::cv_dashboard::volmap_periodic_y } {
          set periodic_string ${periodic_string}"yY"
        }
        if { $::cv_dashboard::volmap_periodic_z } {
          set periodic_string ${periodic_string}"zZ"
        }
        if { [string length ${periodic_string}] > 0 } {
            mol showperiodic $::cv_dashboard::mol $repid ${periodic_string}
        }
        set repname [mol repname $::cv_dashboard::mol $repid]
        lappend repnames $repname
        incr color
      }
    }
    set ::cv_dashboard::volmap_rep($cv) $repnames
  }
}


# Hide volmaps for selected colvars
proc ::cv_dashboard::hide_volmaps_selected {} {
  hide_volmaps [selected_colvars]
}


# Hide volmaps for colvars given as parameters
proc ::cv_dashboard::hide_volmaps { colvars } {
  foreach cv $colvars {
    if { [info exists ::cv_dashboard::volmap_rep($cv)] } {
      foreach r $::cv_dashboard::volmap_rep($cv) {
        mol delrep [mol repindex $::cv_dashboard::mol $r] $::cv_dashboard::mol
      }
      unset ::cv_dashboard::volmap_rep($cv)
    }
  }
}


# Remove all volmap representations
proc ::cv_dashboard::hide_all_volmaps {} {
  foreach { cv repnames } [array get ::cv_dashboard::volmap_rep] {
    foreach r $repnames {
      mol delrep [mol repindex $::cv_dashboard::mol $r] $::cv_dashboard::mol
    }
  }
  array unset ::cv_dashboard::volmap_rep *
}


#################################################################
# Show Colvar gradients
#################################################################


proc ::cv_dashboard::show_gradients { list } {

  foreach cv $list {
    if { ![info exists ::cv_dashboard::grad_objects($cv)] } {
      if { [run_cv colvar $cv set collect_gradient 1] == -1 } { continue }
      run_cv colvar $cv update ;# required to get initial values of gradients
      # Associate empty list of objects to cv to request its update
      set ::cv_dashboard::grad_objects($cv) {}
    }
  }
  update_shown_gradients
}


proc ::cv_dashboard::update_shown_gradients {} {

  set id 0
  set molid $::cv_dashboard::mol
  foreach { cv objs } [array get ::cv_dashboard::grad_objects] {

    # Delete out-of-date graphical objects (arrows)
    foreach obj $objs {
      graphics $molid delete $obj
    }

    # Forget variables that have been deleted (*after* deleting the graphics above)
    if { [lsearch [run_cv list] $cv] == -1 } {
      unset ::cv_dashboard::grad_objects($cv)
      continue
    }

    set atomids [run_cv colvar $cv getatomids]
    if { [llength $atomids] == 0 } {
      # Variable was reinitialized and lost its gradient feature
      if { [run_cv colvar $cv set collect_gradient 1] == -1 } { continue }
      run_cv colvar $cv update
      set atomids [run_cv colvar $cv getatomids]
      # If that didn't work then gradients are not supported
      if { [llength $atomids] == 0 } { continue }
    }

    set grads [run_cv colvar $cv getgradients]
    # Get width if provided in colvar config
    set width 1.
    regexp -nocase -line {^\s*width\s+([\d\.e]*)} [get_cv_config $cv] match width

    set new_objs [draw_vectors $atomids $grads $id $width]
    incr id

    if { [llength $new_objs] > 0 } {
      set ::cv_dashboard::grad_objects($cv) $new_objs
    }
  }
}


proc ::cv_dashboard::hide_gradients {} {
  foreach cv [selected_colvars] {
    if [info exists ::cv_dashboard::grad_objects($cv)] {
      set objs $::cv_dashboard::grad_objects($cv)
      foreach obj $objs { graphics $::cv_dashboard::mol delete $obj }
      run_cv colvar $cv set collect_gradient 0
      unset ::cv_dashboard::grad_objects($cv)
    }
  }
}


proc ::cv_dashboard::hide_all_gradients {} {
  foreach { cv objs } [array get ::cv_dashboard::grad_objects] {
    foreach obj $objs { graphics $::cv_dashboard::mol delete $obj }
    run_cv colvar $cv set collect_gradient 0
  }
  array unset ::cv_dashboard::grad_objects *
}


#################################################################
# Show Bias forces
#################################################################


proc ::cv_dashboard::show_forces { list } {

  foreach bias $list {
    if { ![info exists ::cv_dashboard::force_objects($bias)] } {
      # Associate empty list of objects to cv to request its update
      set ::cv_dashboard::force_objects($bias) {}
    }
  }
  update_shown_forces
}


proc ::cv_dashboard::update_shown_forces {} {

  # Start IDs for force objects after those of gradient objects
  set id [array size ::cv_dashboard::grad_objects]
  set molid $::cv_dashboard::mol
  set atomids [run_cv getatomids]

  foreach { bias objs } [array get ::cv_dashboard::force_objects] {

    # Delete out-of-date graphical objects (arrows)
    foreach obj $objs {
      graphics $molid delete $obj
    }

    # Forget biases that have been deleted (*after* deleting the graphics above)
    if { [lsearch [run_cv list biases] $bias] == -1 } {
      unset ::cv_dashboard::force_objects($bias)
      continue
    }

    # Obtaining forces from each particular biases requires a dedicated Module update
    # enabling forces only from the chosen bias
    foreach b [run_cv list biases] {
      if { $b == $bias } {
        cv bias $b set apply_force 1
      } else {
        cv bias $b set apply_force 0
      }
    }
    # Complete module update is necessary to get atomic bias forces
    run_cv update
    set force_list [run_cv getatomappliedforces]

    set new_objs [draw_vectors $atomids $force_list $id]
    incr id

    if { [llength $new_objs] > 0 } {
      set ::cv_dashboard::force_objects($bias) $new_objs
    }
  }
  # Restore all forces and update Module
  foreach b [run_cv list biases] {
    cv bias $b set apply_force 1
  }
  run_cv update
}

proc ::cv_dashboard::hide_forces {} {
  foreach b [selected_biases] {
    if [info exists ::cv_dashboard::force_objects($b)] {
      set objs $::cv_dashboard::force_objects($b)
      foreach obj $objs { graphics $::cv_dashboard::mol delete $obj }
      unset ::cv_dashboard::force_objects($b)
    }
  }
}


proc ::cv_dashboard::hide_all_forces {} {
  foreach { b objs } [array get ::cv_dashboard::force_objects] {
    foreach obj $objs { graphics $::cv_dashboard::mol delete $obj }
  }
  array unset ::cv_dashboard::force_objects *
}


#################################################################
# Draw vectors as arrows, used by gradient and force displays
#################################################################


proc ::cv_dashboard::draw_vectors { atomids vector_list obj_id { width 1.} } {
  set molid $::cv_dashboard::mol
  set w .cv_dashboard_window

  if { [llength $vector_list] == 0 } { return }

  if  { [llength $vector_list] != [llength $atomids] } {
    puts "Error: mismatched list lengths in ::cv_dashboard::draw_vectors"
    puts $vector_list
    puts $atomids
    return
  }

  # Map gradients list to dictionary
  set maxl2 0.
  foreach vec $vector_list aid $atomids {
    set l2 [veclength2 $vec]
    if { $l2 > $maxl2 } { set maxl2 $l2 }
    set vectors($aid) $vec
  }
  if { $maxl2 < 1e-14 } {
    # Zero vectors, don't even try
    return
  }

  # Now we can safely sort atom IDs
  set atomids [lsort -integer $atomids]
  set sel [atomselect $molid "($::cv_dashboard::sel_text) and (index $atomids)"]
  set coords [$sel get {x y z}]
  # Filter list of atom ids
  set atomids [$sel list]
  $sel delete

  graphics $molid material [.cv_dashboard_window.tabs.settings.material get]

  # Skip colors used by other gradient or force objects
  # avoid colors 0-2
  graphics $molid color [expr $obj_id % 29 + 3 ]

  if { $::cv_dashboard::grad_scale_choice == "scale" } {
    set fact [expr {$::cv_dashboard::grad_scale / $width}]
    set grad_norm [expr {sqrt($maxl2) * $fact}]
    set ::cv_dashboard::grad_norm [round $grad_norm 5]
  } else {
    set fact [expr {$::cv_dashboard::grad_norm / sqrt($maxl2)}]
    set grad_scale [expr {$fact * $width}]
    set ::cv_dashboard::grad_scale [round $grad_scale 5]
  }

  # Create new arrows
  set radius [$w.tabs.settings.grad_radius get]
  set new_objs {}
  foreach start $coords id $atomids {
    set vec [vecscale $fact $vectors($id)]
    set vec_len [veclength $vec]
    # Don't draw zero-length vectors
    if { $vec_len < 1e-2 } { continue }
    set end [vecadd $start $vec]
    if { ${vec_len} > [expr 6.0*${radius}] } {
      # Long arrow : cone length is 3 times radius
      set middle [vecadd $start \
        [vecscale [expr (${vec_len} - 3.0*${radius})/${vec_len}] ${vec}]]
    } else {
      # Short arrow: cap cone length at 1/2 total length
      set middle [vecadd $start [vecscale 0.5 ${vec}]]
    }
    set cyl [graphics $molid cylinder $start $middle radius ${radius} resolution 12]
    set cone [graphics $molid cone $middle $end radius [expr ${radius}*2.0] resolution 12]
    lappend new_objs $cyl $cone
  }
  return $new_objs
}


proc ::cv_dashboard::createVolmapMenu { row } {
  set menu .cv_dashboard_window.tabs.main
  set gridrow $row

  # Volumetric map display settings
  grid [ttk::button $menu.show_volmaps -text "Show volmaps" -command {::cv_dashboard::show_volmaps_selected} -padding "2 0 2 0"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $menu.hide_volmaps -text "Hide volmaps" -command {::cv_dashboard::hide_volmaps_selected} -padding "2 0 2 0"] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $menu.hide_all_volmaps -text "Hide all volmaps" -command {::cv_dashboard::hide_all_volmaps} -padding "2 0 2 0"] -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew

  grid remove $menu.show_volmaps $menu.hide_volmaps $menu.hide_all_volmaps
}


proc ::cv_dashboard::createRotationMenu { row } {

  set menu .cv_dashboard_window.tabs.main

  set gridrow $row
  grid [ttk::button $menu.show_rotation -text "Show rotation" -command { ::cv_dashboard::start_rotation_display {} } -padding "2 0 2 0"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $menu.hide_rotation -text "Hide rotation" -command { ::cv_dashboard::stop_rotation_display } -padding "2 0 2 0"] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew

  grid remove $menu.show_rotation $menu.hide_rotation
}
