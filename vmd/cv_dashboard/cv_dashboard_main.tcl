#################################################################
# Main UI: the dashboard
#################################################################


# Create main window
proc ::cv_dashboard::createWindow {} {

  if {[winfo exists .cv_dashboard_window]} {
    wm deiconify .cv_dashboard_window
    return
  }

  if {[molinfo num] == 0 } {
    tk_messageBox -icon error -title "Colvars Dashboard Error"\
      -message "No molecule loaded. Please load a molecule before running Colvars Dashboard.\n"
    return -1
  }

  set w [toplevel .cv_dashboard_window]
  wm title $w "Colvars dashboard"
  wm protocol $w WM_DELETE_WINDOW {
    # window destructor that removes the trace we put in place, so they don't accumulate
    # if loaded multiple times
    set molid [molinfo top]
    trace remove variable vmd_frame($molid) write ::cv_dashboard::update_frame
    destroy .cv_dashboard_window
  }

  # setup Colvars if not already there
  if [catch { cv version}] {
    run_cv molid top
  }
  set gridrow 0
  grid [ttk::button $w.load -text "Load" -command ::cv_dashboard::load -padding "2 0 2 0"] \
    -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $w.save -text "Save" -command ::cv_dashboard::save -padding "2 0 2 0"] \
    -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $w.reset -text "Reset" -command ::cv_dashboard::reset -padding "2 0 2 0"] \
    -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew

  ttk::treeview $w.cvtable -selectmode extended -show tree
  $w.cvtable configure -column val
  $w.cvtable column #0 -width 50 -stretch 1 -anchor w
  $w.cvtable column val -width 150 -stretch 1 -anchor w
  bind $w.cvtable <e> ::cv_dashboard::edit
  bind $w <Control-a> {
    puts $::cv_dashboard::cvs
    .cv_dashboard_window.cvtable selection set $::cv_dashboard::cvs
  }

  $w.cvtable tag configure parity0 -background white
  $w.cvtable tag configure parity1 -background grey94
  refresh_table

  incr gridrow
  grid $w.cvtable -row $gridrow -column 0 -sticky news -columnspan 3
  grid rowconfigure $w $gridrow -weight 1

  incr gridrow
  grid [ttk::button $w.edit -text "Edit \[e\]" -command ::cv_dashboard::edit -padding "2 0 2 0"] \
    -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $w.add -text "New" -command ::cv_dashboard::add -padding "2 0 2 0"] \
    -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $w.del -text "Delete" -command ::cv_dashboard::del -padding "2 0 2 0"] \
    -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew

  incr gridrow
  grid [ttk::button $w.plot -text "Timeline plot" -command ::cv_dashboard::plot -padding "2 0 2 0"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $w.plot2cv -text "Pairwise plot" -command {::cv_dashboard::plot 2cv} -padding "2 0 2 0"] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $w.refresh -text "Refresh (F5)" -command ::cv_dashboard::refresh_table -padding "2 0 2 0"] -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew

  user add key F5 ::cv_dashboard::refresh_table
  bind $w <F5> ::cv_dashboard::refresh_table

  # Cannot test directly for the presence of scripting methods in the absence of a defined colvar
  # so we test the version number instead
  if {[string compare [run_cv version] "2019-02-07"] >= 0} {
    incr gridrow
    grid [ttk::button $w.show_atoms -text "Show atoms" -command {::cv_dashboard::show_atoms} -padding "2 0 2 0"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
    grid [ttk::button $w.hide_atoms -text "Hide all" -command {::cv_dashboard::hide_atoms} -padding "2 0 2 0"] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  }

  if {[string compare [run_cv version] "2019-03-18"] >= 0} {
    incr gridrow
    grid [ttk::button $w.show_gradients -text "Show gradients" -command {::cv_dashboard::show_gradients [::cv_dashboard::selected_colvars]} \
      -padding "2 0 2 0"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
    grid [ttk::button $w.hide_gradients -text "Hide all" -command {::cv_dashboard::hide_gradients} -padding "2 0 2 0"] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  }

  incr gridrow
  grid [label $w.frameTxt -text "Frame:"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [label $w.frame -textvariable ::cv_dashboard::current_frame] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [ttk::checkbutton $w.trackFrame -text "Track VMD frame" -command ::cv_dashboard::change_track_frame -variable ::cv_dashboard::track_frame] \
    -row $gridrow -column 2  -pady 2 -padx 2 -sticky nsew
  change_track_frame ;# activate tracking if necessary

  grid columnconfigure $w 0 -weight 1
  grid columnconfigure $w 1 -weight 1
  grid columnconfigure $w 2 -weight 1

  return $w
}


# Refresh the table with a list of existing CVs and their values
proc ::cv_dashboard::refresh_table {} {

  set w .cv_dashboard_window
  set ::cv_dashboard::cvs [run_cv list]

  foreach i [$w.cvtable children {}] {
    $w.cvtable delete $i
  }

  # Get fresh coordinates from VMD
  run_cv update

  set parity 1
  foreach c $::cv_dashboard::cvs {
    $w.cvtable insert {} end -id $c -text $c
    $w.cvtable tag add cvname [list $c]
    # tag odd and even rows for alternating background colors
    set parity [expr 1-$parity]
    $w.cvtable tag add "parity$parity" [list $c]

    set val [run_cv colvar $c update]
    $w.cvtable set $c val [format_value $val]

    set size [llength $val]
    if { $size > 1 } {
      for {set i 1} {$i <= $size} {incr i} {
        set n "${c} ${i}"
        $w.cvtable insert $c end -id $n -text $n
        # all colvar/scalar comp items are tagged with tag cvname, to define common bindings on them
        $w.cvtable tag add cvname [list $n]
        $w.cvtable tag add "parity$parity" [list $n]

        $w.cvtable set $n val [format_value [lindex $val [expr $i-1]]]
      }
    }
  }

  # Remove deleted variables from gradient display
  set tmp {}
  foreach cv $::cv_dashboard::grad_cvs {
    if { [lsearch $::cv_dashboard::cvs $cv] > -1 } { lappend tmp $cv }
  }
  set ::cv_dashboard::grad_cvs $tmp
  show_gradients $::cv_dashboard::grad_cvs

  update_frame internal [molinfo top] w
}


# Refresh the table with a list of colvar values
proc ::cv_dashboard::refresh_values {} {
  run_cv update
  set w .cv_dashboard_window

  foreach c [$w.cvtable children {}] {
    set val [run_cv colvar $c update]
    $w.cvtable set $c val [format_value $val]

    set size [llength $val]
    if { $size > 1 } {
      for {set i 1} {$i <= $size} {incr i} {
        set n "${c} ${i}"
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
      # reduce salar components to their parent CV
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
  global vmd_frame

  set molid [molinfo top]
  if {$::cv_dashboard::track_frame} {
    trace add variable vmd_frame($molid) write ::cv_dashboard::update_frame
    update_frame internal [molinfo top] w
  } else {
    trace remove variable vmd_frame($molid) write ::cv_dashboard::update_frame
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
    run_cv configfile $path
    refresh_table
    # Save directory for next invocation of this dialog
    set ::cv_dashboard::config_dir [file dirname $path]
  }
}


# Save config of colvars to file (can we do it w/ biases ? need bias type keyword)
proc ::cv_dashboard::save {} {
  if { [info exists ::cv_dashboard::config_dir] } {
    set path [tk_getSaveFile -filetypes {{"Colvars cfg" .in} {"Colvars cfg" .colvars} {"Gromacs Colvars cfg" .dat} {"All files" *}} \
        -initialdir $::cv_dashboard::config_dir]
  } else {
    set path [tk_getSaveFile -filetypes {{"Colvars cfg" .in} {"Colvars cfg" .colvars} {"Gromacs Colvars cfg" .dat} {"All files" *}} \
        -initialdir [pwd]]
  }

  if [string compare $path ""] {
    # Save directory for next invocation of this dialog
    set ::cv_dashboard::config_dir [file dirname $path]
    set cfg ""
    foreach c [run_cv list] {
        append cfg "colvar {"
        append cfg [run_cv colvar $c getconfig]
        append cfg "}\n\n"
    }

    set o [open $path w]
    puts $o $cfg
    close $o
  }
}


# Delete currently selected colvars
proc ::cv_dashboard::del {} {
  foreach c [selected_colvars] {
    # workaround bug in colvars pre-2018-07-02
    if {[string compare [run_cv version] "2018-07-02"] == -1} {
      foreach b [run_cv list biases] {
        if [string match *colvars*$c* [run_cv bias $b getconfig]] {
          puts "Bug workaround: deleting bias $b"
          run_cv bias $b delete
        }
      }
    }
    run_cv colvar $c delete
  }
  refresh_table
}


# Reset cvm: hard delete allows for reconstructing the module after changing top molecule.
proc ::cv_dashboard::reset {} {
  run_cv delete
  run_cv molid top
  refresh_table
}


#################################################################
# Display atoms
#################################################################


# Display atoms in groups for selected colvars
proc ::cv_dashboard::show_atoms {} {
  set color 0
  set ci 0
  foreach c [selected_colvars] {
    incr ci
    set all_groups [run_cv colvar $c getatomgroups]
    # Remove characters such as <> which are parsed as special in VMD selection texts
    set sanitized_cvname [regsub -all {[^a-zA-Z0-9_@]} $c {} ]
    set i 0
    foreach list $all_groups {
      incr i
      # dummyAtoms will return empty lists
      if {[llength $list] > 0} {
        set group "${sanitized_cvname}_group_${i}"
        # resolve ambiguous names due to colvar name sanitization
        while {[lsearch $::cv_dashboard::macros $group] > -1} {
          append sanitized_cvname "_"
          set group "${sanitized_cvname}_group_${i}"
        }
        atomselect macro $group "index $list"
        lappend ::cv_dashboard::macros $group
        mol color ColorID $color
        mol representation VDW 0.5 12.
        mol selection "$group"
        mol material Opaque
        mol addrep top
        set repid [expr [molinfo top get numreps] - 1]
        lappend ::cv_dashboard::repnames [mol repname top $repid]
        incr color
      }
    }
  }
}


# Remove atom representations
proc ::cv_dashboard::hide_atoms {} {
  foreach r $::cv_dashboard::repnames {
    mol delrep [mol repindex top $r] top
  }
  set ::cv_dashboard::repnames {}
  set ::cv_dashboard::macros [lsort -unique $::cv_dashboard::macros]
  foreach m $::cv_dashboard::macros {
    atomselect delmacro $m
  }
  set ::cv_dashboard::macros {}
}


#################################################################
# Show Colvar gradients
#################################################################


proc ::cv_dashboard::show_gradients { list } {

  foreach cv $list {
    if { [run_cv colvar $cv set "collect gradient" 1] == -1 } { continue }
    run_cv colvar $cv update ;# required to get inital values of gradients
    if { [lsearch $::cv_dashboard::grad_cvs $cv] == -1 } { lappend ::cv_dashboard::grad_cvs $cv }
  }
  update_shown_gradients
}


proc ::cv_dashboard::update_shown_gradients {} {

  foreach i $::cv_dashboard::grad_objects {
    graphics top delete $i
  }
  set colorid 3 ;# avoid very common or less visible colors blue, red, gray
  foreach cv $::cv_dashboard::grad_cvs {
    set atomids [run_cv colvar $cv getatomids]
    if { [llength $atomids] == 0 } { continue }
    set grads [run_cv colvar $cv getgradients]
    set sel [atomselect top "index $atomids"]
    set coords [$sel get {x y z}]
    $sel delete

    set desired_max_length 5.
    graphics top color [expr $colorid % 32]
    incr colorid

    set maxl2 0.

    foreach g $grads {
      set l2 [veclength2 $g]
      if { $l2 > $maxl2 } { set maxl2 $l2 }
    }
    if {$maxl2 < 1e-10} { continue }
    set fact [expr {$desired_max_length / sqrt($maxl2)}]

    foreach i $atomids start $coords g $grads {
      set vec [vecscale $fact $g]
      set end [vecadd $start $vec]
      set middle [vecadd $start [vecscale 0.9 $vec]]
      set cyl [graphics top cylinder $start $middle radius 0.15]
      set cone [graphics top cone $middle $end radius 0.25]
      lappend ::cv_dashboard::grad_objects $cyl $cone
    }
  }
}


proc ::cv_dashboard::hide_gradients {} {

  foreach i $::cv_dashboard::grad_objects {
    graphics top delete $i
  }
  set ::cv_dashboard::grad_objects {}
  foreach cv $::cv_dashboard::grad_cvs {
    run_cv colvar $cv set "collect gradient" 0
  }
  set ::cv_dashboard::grad_cvs {}
}
