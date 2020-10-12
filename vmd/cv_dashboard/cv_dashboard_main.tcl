#################################################################
# Main UI: the dashboard
#################################################################

# -*- tcl-indent-level: 2; -*-


# Create main window
proc ::cv_dashboard::createWindow {} {

  set w [toplevel .cv_dashboard_window]
  wm title $w "Colvars dashboard"
  wm protocol $w WM_DELETE_WINDOW { ::cv_dashboard::quit }

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
  grid [ttk::button $w.save -text "Save colvars" -command ::cv_dashboard::save -padding "2 0 2 0"] \
    -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $w.reset -text "Reset" -command ::cv_dashboard::reset -padding "2 0 2 0"] \
    -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew

  ttk::treeview $w.cvtable -selectmode extended -show tree
  $w.cvtable configure -column val
  $w.cvtable column #0 -width 50 -stretch 1 -anchor w
  $w.cvtable column val -width 150 -stretch 1 -anchor w

  bind $w.cvtable <Button-3> {::cv_dashboard::cvContextMenu %x %y %X %Y}
  bind $w.cvtable <Button-1> {::cv_dashboard::cvTableClicked %x %y}

  bind $w <Control-e> ::cv_dashboard::edit
  bind $w <Control-a> { .cv_dashboard_window.cvtable selection set $::cv_dashboard::cvs }
  bind $w <Control-n> ::cv_dashboard::add

  event add <<keyb_enter>> <Return>   ;# Combine Return and keypad-Enter into a single virtual event
  event add <<keyb_enter>> <KP_Enter>

  if { [info patchlevel] != "8.5.6" } {
    $w.cvtable tag configure parity0 -background white
    $w.cvtable tag configure parity1 -background grey94
  }
  refresh_table

  incr gridrow
  grid $w.cvtable -row $gridrow -column 0 -sticky news -columnspan 3
  grid rowconfigure $w $gridrow -weight 1

  # Editing
  incr gridrow
  grid [ttk::button $w.edit -text "Edit \[Ctrl-e\]" -command ::cv_dashboard::edit -padding "2 0 2 0"] \
    -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $w.add -text "New \[Ctrl-n\]" -command ::cv_dashboard::add -padding "2 0 2 0"] \
    -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $w.del -text "Delete" -command ::cv_dashboard::del -padding "2 0 2 0"] \
    -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew

  # Plots
  incr gridrow
  grid [ttk::button $w.plot -text "Timeline plot" -command ::cv_dashboard::plot -padding "2 0 2 0"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $w.plot2cv -text "Pairwise plot" -command {::cv_dashboard::plot 2cv} -padding "2 0 2 0"] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $w.refresh -text "Refresh \[F5\]" -command ::cv_dashboard::refresh_table -padding "2 0 2 0"] -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew

  user add key F5 ::cv_dashboard::refresh_table
  bind $w <F5> ::cv_dashboard::refresh_table

  # Add trajectory animation bindings to Dashboard and VMD's OpenGL window
  traj_animation_bindings $w
  user add key Left  { ::cv_dashboard::chg_frame -1 }
  user add key Right { ::cv_dashboard::chg_frame 1 }
  user add key Home  { ::cv_dashboard::chg_frame start }
  user add key End   { ::cv_dashboard::chg_frame end }

  incr gridrow
  grid [ttk::separator $w.sep0 -orient horizontal] -row $gridrow -column 0 -columnspan 3 -pady 5 -sticky ew

  # Atom group display
  incr gridrow
  grid [ttk::button $w.show_atoms -text "Show atoms" -command {::cv_dashboard::show_atoms_selected} -padding "2 0 2 0"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $w.hide_atoms -text "Hide atoms" -command {::cv_dashboard::hide_atoms_selected} -padding "2 0 2 0"] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $w.hide_all_atoms -text "Hide all atoms" -command {::cv_dashboard::hide_all_atoms} -padding "2 0 2 0"] -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew

  # Gradient display
  incr gridrow
  grid [ttk::button $w.show_gradients -text "Show gradients" -command {::cv_dashboard::show_gradients [::cv_dashboard::selected_colvars]} \
    -padding "2 0 2 0"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $w.hide_gradients -text "Hide gradients" -command {::cv_dashboard::hide_gradients} -padding "2 0 2 0"] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $w.hide_all_gradients -text "Hide all grads" -command {::cv_dashboard::hide_all_gradients} -padding "2 0 2 0"] -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew

  incr gridrow
  grid [ttk::separator $w.sep2 -orient horizontal] -row $gridrow -column 0 -columnspan 3 -pady 5 -sticky ew

  # Toggle advanced representations

  incr gridrow
  grid [ttk::button $w.volmapBtn -text "Volmap menu" -command ::cv_dashboard::toggleVolmapMenu -padding "2 0 2 0"] -row $gridrow -column 0 -columnspan 1 -pady 2 -padx 2 -sticky nsew

  incr gridrow
  grid [ttk::separator $w.sep3 -orient horizontal] -row $gridrow -column 0 -columnspan 3 -pady 5 -sticky ew

  # Molecule
  incr gridrow
  grid [label $w.molTxt -text "Molecule:"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  ttk::combobox $w.mol -justify left -state readonly
  $w.mol configure -values [molinfo list]
  if { $::cv_dashboard::mol != -1 } {
    $w.mol set $::cv_dashboard::mol
  }
  trace add variable ::vmd_initialize_structure write ::cv_dashboard::update_mol_list
  grid $w.mol -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  bind $w.mol <<ComboboxSelected>> ::cv_dashboard::change_mol

  # Frame display and track checkbox
  incr gridrow
  grid [label $w.frameTxt -text "Frame:"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [label $w.frame -textvariable ::cv_dashboard::current_frame] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [ttk::checkbutton $w.trackFrame -text "Track VMD frame" -command ::cv_dashboard::change_track_frame -variable ::cv_dashboard::track_frame] \
    -row $gridrow -column 2  -pady 2 -padx 2 -sticky nsew
  change_track_frame ;# activate tracking if necessary

  incr gridrow
  grid [ttk::button $w.settingsBtn -text "Show advanced settings" -command ::cv_dashboard::toggleSettings -padding "2 0 2 0"] -row $gridrow -column 0 -columnspan 3 -pady 2 -padx 2 -sticky nsew

  # Create and hide volmap menu
  createVolmapMenu

  # Create and hide Settings window to create all associated variables
  createSettingsWindow

  grid columnconfigure $w 0 -weight 1
  grid columnconfigure $w 1 -weight 1
  grid columnconfigure $w 2 -weight 1

  return $w
}


# Open or close the volmap sub-panel
proc ::cv_dashboard::toggleVolmapMenu {} {
  set w .cv_dashboard_window
  if { $::cv_dashboard::volmap_menu_shown } {
    set ::cv_dashboard::volmap_menu_shown false
    grid remove $w.volmap_menu
  } else {
    set ::cv_dashboard::volmap_menu_shown true
    grid $w.volmap_menu
  }
}


# Open or close the Settings panel
proc ::cv_dashboard::toggleSettings {} {
  set w .cv_dashboard_window

  if { $::cv_dashboard::settings_shown } {
    set ::cv_dashboard::settings_shown false
    $w.settingsBtn configure -text "Show advanced settings"
    grid remove $w.settings
  } else {
    set ::cv_dashboard::settings_shown true
    $w.settingsBtn configure -text "Hide advanced settings"
    grid $w.settings
  }
}


# Display context menu about specific colvar
# Takes coordinates within widget and within window
proc ::cv_dashboard::cvContextMenu { x y wX wY } {
  set w .cv_dashboard_window
  set menu $w.cvMenu

  set cvs [selected_colvars]
  # Add any colvar under mouse cursor
  foreach cv [$w.cvtable identify item $x $y] {
    if { [lsearch $cvs $cv] == -1 } {
      # reduce scalar components to their parent vector CV
      lappend cvs [lindex $cv 0]
    }
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
    $menu add command -label New -command ::cv_dashboard::add
  } else {
    if { [llength $volmaps] > 0 } {
      $menu add command -label "Show volmap" -command [list ::cv_dashboard::show_volmaps $volmaps]
      $menu add command -label "Hide volmap" -command [list ::cv_dashboard::hide_volmaps $volmaps]
    }
    if { [llength $rotations] > 0 } {
      $menu add command -label "Show rotation" -command [list ::cv_dashboard::start_rotation_display $rotations]
      $menu add command -label "Hide rotation" -command [list ::cv_dashboard::stop_rotation_display]
    }
    $menu add command -label Edit -command [list ::cv_dashboard::edit false $cvs]
    $menu add command -label Delete -command [list ::cv_dashboard::del $cvs]
  }
  tk_popup $menu $wX $wY
}


# Takes coordinates within widget
proc ::cv_dashboard::cvTableClicked { x y } {
  set w .cv_dashboard_window
  set menu $w.cvMenu

  # colvar under mouse
  set cv [$w.cvtable identify item $x $y]

  # Deselect all if clicked on nothing
  if { [llength $cv] == 0 } {
    $w.cvtable selection set [list]
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

G. Fiorin, M. L. Klein, and J. Hénin. Using collective variables to drive molecular dynamics simulations. Mol. Phys., 111(22-23):3345–3362, 2013
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
}


# Format colvar values returned by cv for display in table
proc ::cv_dashboard::format_value val {
  if {[llength $val] == 1} {
    return [format "%.4g" $val]
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
    set path [tk_getOpenFile -filetypes {{"Colvars cfg" .in} {"Colvars cfg" .colvars} {"Gromacs Colvars cfg" .dat} {"All files" *}} \
        -initialdir $::cv_dashboard::config_dir]
  } else {
    set path [tk_getOpenFile -filetypes {{"Colvars cfg" .in} {"Colvars cfg" .colvars} {"Gromacs Colvars cfg" .dat} {"All files" *}} \
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
    set path [tk_getSaveFile -confirmoverwrite false -filetypes {{"Colvars cfg" .in} {"Colvars cfg" .colvars} {"Gromacs Colvars cfg" .dat} {"All files" *}} \
        -initialdir $::cv_dashboard::config_dir]
  } else {
    set path [tk_getSaveFile -confirmoverwrite false -filetypes {{"Colvars cfg" .in} {"Colvars cfg" .colvars} {"Gromacs Colvars cfg" .dat} {"All files" *}} \
        -initialdir [pwd]]
  }

  if [string compare $path ""] {
    # Save directory for next invocation of this dialog
    set ::cv_dashboard::config_dir [file dirname $path]
    if {$::cv_dashboard::units == ""} {
      set cfg ""
    } else {
      set cfg "units $::cv_dashboard::units\n\n"
    }
    set indexFiles [list]
    catch { set indexFiles [cv listindexfiles] }
    foreach ndx $indexFiles {
      append cfg "indexFile $ndx\n"
    }
    foreach cv [run_cv list] {
        append cfg "colvar {[get_config $cv]}\n\n"
    }
    if [file exists $path] {
      set answer [tk_messageBox -icon warning -type okcancel -title "Warning: file overwrite"\
        -message "Saving to [file tail $path] will overwrite it.\n\nOnly definitions of collective variables - not biases - will be saved!\n\nProceed?"]
      if { $answer == "cancel" } { return }
    }
    set o [open $path w]
    puts $o $cfg
    close $o
  }
}


# Delete currently selected colvars
proc ::cv_dashboard::del { {cvs "" } } {
  if { [llength $cvs] < 1 } {
    set cvs [selected_colvars]
  }
  foreach cv $cvs {
    # workaround bug in colvars pre-2018-07-02
    if {[string compare [run_cv version] "2018-07-02"] == -1} {
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
  }

  # Run cv delete silently to be less verbose if module was already deleted
  catch { cv delete }
  if { $molid != -1 } {
    run_cv molid $molid
  }
  set ::cv_dashboard::colvar_configs [dict create]
  refresh_table
  refresh_units
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
        atomselect macro $group "index $list"
        mol color ColorID $color
        mol representation VDW 0.5 12.
        mol selection "$group"
        mol material Opaque
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
    set cv_volmaps [run_cv colvar $cv getvolmapids]
    set repnames {}
    foreach volid $cv_volmaps {
      if { $volid >= 0 } {
        set threshold [expr $::cv_dashboard::volmap_contour * [vecsum \
          [voltool info minmax -mol $::cv_dashboard::mol -vol ${volid}]]]
        mol color ColorID $color
        mol selection all  ;# Must provide some selection text
        mol representation Isosurface ${threshold} ${volid} 2 0 0 1
        mol material [$w.volmap_menu.volmap_material get]
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

  set colorid 3 ;# avoid very common or less visible colors blue, red, gray
  set molid $::cv_dashboard::mol
  foreach { cv objs } [array get ::cv_dashboard::grad_objects] {

    # Delete out-of-date graphical objects (arrows)
    foreach obj $objs {
      graphics $molid delete $obj
    }

    # Forget variables that have been deleted (*after* deleting the graphics above)
    if { [lsearch $::cv_dashboard::cvs $cv] == -1 } {
      unset ::cv_dashboard::grad_objects($cv)
      continue
    }

    set atomids [run_cv colvar $cv getatomids]
    if { [llength $atomids] == 0 } {
      # Variable was reinitialized and lost its gradient feature
      if { [run_cv colvar $cv set collect_gradient 1] == -1 } { continue }
      # If that didn't work then gradients are not supported
      if { [llength $atomids] == 0 } { continue }
      run_cv colvar $cv update
    }

    set grads [run_cv colvar $cv getgradients]
    if { [llength $grads] == 0 } { continue }

    set sel [atomselect $molid "index $atomids"]
    set coords [$sel get {x y z}]
    $sel delete

    # Loop through colorids (only in this run of the proc though)
    graphics $molid color [expr $colorid % 32]
    incr colorid

    set maxl2 0.
    foreach g $grads {
      set l2 [veclength2 $g]
      if { $l2 > $maxl2 } { set maxl2 $l2 }
    }
    if {$maxl2 < 1e-14} { continue } ;# Zero gradient, don't even try

    # Get width if provided in colvar config
    set width 1.
    regexp -nocase -lineanchor {^\s*width\s([\d\.e]*)} [get_config $cv] match width

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
    set new_objs {}
    foreach start $coords g $grads {
      set vec [vecscale $fact $g]
      set end [vecadd $start $vec]
      set middle [vecadd $start [vecscale 0.9 $vec]]
      set cyl [graphics $molid cylinder $start $middle radius 0.1 resolution 12]
      set cone [graphics $molid cone $middle $end radius 0.2 resolution 12]
      lappend new_objs $cyl $cone
    }
    set ::cv_dashboard::grad_objects($cv) $new_objs
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


proc ::cv_dashboard::createVolmapMenu { } {

  set w .cv_dashboard_window

  set menu $w.volmap_menu
  grid [frame $menu] -column 0 -columnspan 3 -sticky nsew

  set gridrow 0

  # Volumetric map display settings
  incr gridrow
  grid [ttk::button $menu.show_volmaps -text "Show volmaps" -command {::cv_dashboard::show_volmaps_selected} -padding "2 0 2 0"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $menu.hide_volmaps -text "Hide volmaps" -command {::cv_dashboard::hide_volmaps_selected} -padding "2 0 2 0"] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $menu.hide_all_volmaps -text "Hide all volmaps" -command {::cv_dashboard::hide_all_volmaps} -padding "2 0 2 0"] -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew

  incr gridrow
  grid [label $menu.volmap_material_text -text "Volmap material:"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  ttk::combobox $menu.volmap_material -justify left -state readonly
  $menu.volmap_material configure -values [material list]
  grid $menu.volmap_material -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  $menu.volmap_material set "Opaque"

  incr gridrow
  grid [label $menu.volmap_contour_text -text "Contour level:"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [tk::entry $menu.volmap_contour -textvariable ::cv_dashboard::volmap_contour] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  set ::cv_dashboard::volmap_contour 0.5
  grid [label $menu.volmap_contour_unit -text "(% of min-max range)"] -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew

  incr gridrow
  grid [ttk::checkbutton $menu.volmap_periodic_x -text "+/-X images" -variable ::cv_dashboard::volmap_periodic_x] \
    -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::checkbutton $menu.volmap_periodic_y -text "+/-Y images" -variable ::cv_dashboard::volmap_periodic_y] \
    -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [ttk::checkbutton $menu.volmap_periodic_z -text "+/-Z images" -variable ::cv_dashboard::volmap_periodic_z] \
    -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew
  set ::cv_dashboard::volmap_periodic_x 0
  set ::cv_dashboard::volmap_periodic_y 0
  set ::cv_dashboard::volmap_periodic_z 0

  incr gridrow
  grid [ttk::separator $w.sepvolmap -orient horizontal] -row $gridrow -column 0 -columnspan 3 -pady 5 -sticky ew

  grid columnconfigure $menu 0 -weight 1
  grid columnconfigure $menu 1 -weight 1
  grid columnconfigure $menu 2 -weight 1

  grid remove $menu
  set ::cv_dashboard::volmap_menu_shown false
}
