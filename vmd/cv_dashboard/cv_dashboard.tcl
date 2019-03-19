# Colvars Dashboard -- based on the Colvars Module for VMD
# Jérôme Hénin <henin@ibpc.fr> 2018

# Usage:
# source cv_dashboard.tcl
# cv_dashboard (to reopen window)

# Design principles:
# - take advantage of colvars/VMD binding for maximum user interaction
# - hide the colvars config text from user, instead expose colvar, names and values
# - do not try to parse the colvars config (let the Colvars Module do it)
#   to avoid coming up with an incompatible parser

# This plugin only acts on the "top" molecule
# which is most consistent for trajectory animation (determined by the frame number of top mol)

# TODO PRIORITY:
# - currently only way to handle harmonic walls is legacy, bc separate biases are not accessible!
# - documentation -> link to section in HTML VMD/Colvars on github.io

# TODO:
# - properly calculate position of cursor in plot when not all the plot is visible (resized window)
# - retain comments in config strings (needs upstream support in colvars)
# - Retain the whole config string minus colvars, copy to output? Problem: not visible in interface
# - histograms - either directly from plot window

# TODO Multiplot:
# - handle several windows at once - at least one pairwise, one timeline
# - integrate interactive hacks into interface
# - display pairwise traj on top of known 2D data (eg. FE surface)
# - embed inside main window?

# TODO maybe:
# - index group builder
# - graphical representations such as rotation_display


package provide cv_dashboard 1.0

namespace eval ::cv_dashboard {
  # General UI state
  variable current_frame 0  ;# linked to frame display
  variable cvs {}           ;# linked to colvar table
  variable track_frame 1    ;# start tracking by default

  # State variables for config editor
  variable being_edited
  variable backup_cfg
  variable filetype "atomsFile"

  # Handle to keep track of interactive plot
  variable plothandle
  variable plottype     ;# either timeline or 2cv

  variable repnames {}      ;# representations created by us
  variable macros {}        ;# macros created by us
  variable grad_objects {}  ;# ids of graphical objects displaying gradients
  variable grad_cvs {}      ;# cvs whose gradients are shown

  variable template_dir
  variable template_base_dir
  # Use template dir if full distribution is provided and path is known
  if [info exists ::env(CV_DASHBOARD_DIR)] {
    set template_base_dir ${::env(CV_DASHBOARD_DIR)}/templates
    set template_dir $template_base_dir
  } else {
    set template_dir [pwd]
  }
}


#################################################################
# Main UI: the dashboard
#################################################################


proc cv_dashboard {} {
  return [eval ::cv_dashboard::createWindow]
}


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
    wm destroy .cv_dashboard_window
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
  bind $w <Control-a> { .cv_dashboard_window.cvtable selection set $::cv_dashboard::cvs }

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
  grid [ttk::button $w.refresh -text "Refresh table" -command ::cv_dashboard::refresh_table -padding "2 0 2 0"] -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew

  # Cannot test directly for the presence of getatomids method in the absence of a defined colvar
  # so we test the version number instead
  if {[string compare [run_cv version] "2019-02-07"] >= 0} {
    incr gridrow
    grid [ttk::button $w.show_atoms -text "Show atoms" -command {::cv_dashboard::show_atoms} -padding "2 0 2 0"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
    grid [ttk::button $w.hide_atoms -text "Hide all" -command {::cv_dashboard::hide_atoms} -padding "2 0 2 0"] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  }

  if {[string compare [run_cv version] "2019-03-18"] >= 0} {
    incr gridrow
    grid [ttk::button $w.show_gradients -text "Show gradients" -command {::cv_dashboard::show_gradients} -padding "2 0 2 0"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
    grid [ttk::button $w.hide_gradients -text "Hide all" -command {::cv_dashboard::hide_gradients} -padding "2 0 2 0"] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
    grid [label $w.gradfactor -textvariable ::cv_dashboard::grad_factor] -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew
    # TODO add labek and text frame within 3rd column (sub frame?)
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
    set all_groups [run_cv colvar $c getatomids]
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
        mol representation VDW 1.000000 12.000000
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


proc draw_arrow {mol start end} {
  set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]]
  set ids [list [graphics $mol cylinder $start $middle radius 0.15]]
  lappend ids [graphics $mol cone $middle $end radius 0.25]
  return $ids
}


proc ::cv_dashboard::show_gradients {} {

  foreach cv [selected_colvars] {
    if { [lsearch $::cv_dashboard::grad_cvs $cv] > -1 } { continue }
    if { [run_cv colvar $cv set "collect gradient" 1] == -1 } { continue }
    run_cv colvar $cv update ;# required to get inital values of gradients
    lappend ::cv_dashboard::grad_cvs $cv
  }
  update_shown_gradients
}


proc ::cv_dashboard::update_shown_gradients {} {

  foreach i $::cv_dashboard::grad_objects {
    graphics top delete $i
  }
  set colorid 1
  foreach cv $::cv_dashboard::grad_cvs {

    # For atom groups, could detect identical gradients on atoms, and display based on COM
    set atomids [run_cv colvar $cv getatomids_flat]
    set grads [run_cv colvar $cv getgradients]
    set sel [atomselect top "index $atomids"]
    set coords [$sel get {x y z}]
    $sel delete

    set desired_max_length 5.
    # set colors { yellow red white purple }
    # graphics top color [lindex $colors [expr $colorid % [llength $colors]]]
    graphics top color [expr $colorid % 32]
    incr colorid

    set maxl2 0.

    foreach g $grads {
      set l2 [veclength2 $g]
      if { $l2 > $maxl2 } { set maxl2 $l2 }
    }
    if {$maxl2 < 1e-10} { continue }
    set fact [expr {$desired_max_length / sqrt($maxl2)}]

    foreach i $atomids r $coords g $grads {
      set end [vecadd $r [vecscale $fact $g]]
      set ids [draw_arrow top $r $end]
      lappend ::cv_dashboard::grad_objects {*}$ids
    }
  }
}


proc ::cv_dashboard::hide_gradients {} {

  foreach i $::cv_dashboard::grad_objects {
    graphics top delete $i
  }
  set ::cv_dashboard::grad_objects {}
  set ::cv_dashboard::grad_cvs {}
}


#################################################################
# General utilities
#################################################################


# Call the "cv" interface to Colvars, catching errors and displaying them to the user
proc run_cv args  {
  if [ catch { cv {*}$args } res ] {
    tk_messageBox -icon error -title "Colvars error" -parent .cv_dashboard_window\
      -message "Error running command:\n$args" -detail "$res\n\nSee console for further details."
    return -1
  }
  return $res
}


# Callback to update CV Dashboard when VMD's top molecule changes to new frame
proc ::cv_dashboard::update_frame { name molid op } {
  # name == vmd_frame
  # molid == molecule id of the newly changed frame
  # op == w
  global vmd_frame

  if { $molid != [molinfo top] } {
    return
  }
  set f [molinfo $molid get frame]
  set ::cv_dashboard::current_frame $f

  # set Colvars Module to requested frame
  run_cv frame $f
  # refresh dashboard table
  refresh_values
  # refresh the frame marker in the plot
  display_marker $f
  # refresh displayed CV gradients
  update_shown_gradients
}


#################################################################
# Editor window
#################################################################


# Edit new colvar config
proc ::cv_dashboard::add {} {
  edit true
}


# Colvar config editor window
proc ::cv_dashboard::edit { {add false} } {
  set cfg ""

  if $add {
    # do not remove existing vars
    set cvs {}
    if { [info exists ::cv_dashboard::template_base_dir] } {
      # Open "official" colvar template
      set in [open ${::cv_dashboard::template_base_dir}/colvar.in r]
      set cfg [read $in]
      close $in
    } else {
      # Provide simple template
      set cfg "# You can edit or replace the example colvar config below.\n\
  colvar {\n  name d\n  distance {\n    group1 { atomNumbers 1 2 }\n    group2 { atomNumbers 3 4 }\n  }\n}\n"
    }
  } else {
    set cvs [selected_colvars]
    if {[llength $cvs] == 0} {
      # If no selection, edit all variables
      set cvs $::cv_dashboard::cvs
    }
    foreach c $cvs {
      append cfg "colvar {"
      append cfg [run_cv colvar $c getconfig]
      append cfg "}\n\n"
    }
    set ::cv_dashboard::backup_cfg $cfg
  }
  set w .cv_dashboard_window

  catch { destroy $w.editor }
  set editor [toplevel $w.editor]
  wm title $editor "Colvar config editor"

  # Left frame: utility buttons
  frame $w.editor.fl
  set gridrow 0

  labelframe  $w.editor.fl.docs -bd 2 -text "Online documentation" -padx 2 -pady 2
  set docs $w.editor.fl.docs
  ttk::button $docs.onlinedoc1 -text "Collective variables" -padding "4 2 4 2"\
    -command [list ::cv_dashboard::invokeBrowser "http://colvars.github.io/colvars-refman-vmd/colvars-refman-vmd.html#sec:colvar"]
  ttk::button $docs.onlinedoc3 -text "Basis functions (components)" -padding "4 2 4 2"\
    -command [list ::cv_dashboard::invokeBrowser "http://colvars.github.io/colvars-refman-vmd/colvars-refman-vmd.html#sec:cvc_list"]
  ttk::button $docs.onlinedoc2 -text "Atom groups" -padding "4 2 4 2"\
    -command [list ::cv_dashboard::invokeBrowser "http://colvars.github.io/colvars-refman-vmd/colvars-refman-vmd.html#sec:colvar_atom_groups"]

  grid $docs.onlinedoc1 -row $gridrow -column 0 -pady 5 -padx 2 -sticky nsew
  grid $docs.onlinedoc2 -row $gridrow -column 1 -pady 5 -padx 2 -sticky nsew
  grid $docs.onlinedoc3 -row $gridrow -column 2 -pady 5 -padx 2 -sticky nsew
  grid columnconfigure $docs 0 -weight 1
  grid columnconfigure $docs 1 -weight 1
  grid columnconfigure $docs 2 -weight 1

  labelframe  $w.editor.fl.helpers -bd 2 -text "Editing helpers" -padx 2 -pady 2
  set helpers $w.editor.fl.helpers
  ############# Templates #########################################
  tk::label $helpers.template_label -text "Insert template:"
  ttk::button $helpers.insert_template -text "Pick template file" \
    -command [list ::cv_dashboard::insert_template] -padding "2 0 2 0"

  grid $helpers.template_label -row $gridrow -column 0 -pady 2 -padx 2
  grid $helpers.insert_template -row $gridrow -column 1 -columnspan 2 -sticky ew -pady 2 -padx 2
  incr gridrow

  ############# Atoms from seltext ################################
  tk::label $helpers.seltext_label -text "Atoms from selection text:"
  tk::entry $helpers.seltext -bg white
  # Bind Return key in seltext entry to proc creating the atomNumbers line
  bind $helpers.seltext <Return> "::cv_dashboard::atoms_from_sel textbox"
  ttk::button $helpers.fromsel -text "Insert \[Enter\]" \
    -command "::cv_dashboard::atoms_from_sel textbox" -padding "2 0 2 0"

  grid $helpers.seltext_label -row $gridrow -column 0 -pady 2 -padx 2
  grid $helpers.seltext -row $gridrow -column 1 -sticky ew -pady 2 -padx 2
  grid $helpers.fromsel -row $gridrow -column 2 -pady 2 -padx 2
  incr gridrow

  ############# Atoms from representation ################################
  tk::label $helpers.rep_label -text "Atoms from representation:"
  ttk::combobox $helpers.reps -justify left -state readonly
  ttk::button $helpers.refresh_reps -text "Refresh list" -command ::cv_dashboard::refresh_reps
  bind $helpers.reps <<ComboboxSelected>> "::cv_dashboard::atoms_from_sel reps"

  grid $helpers.rep_label -row $gridrow -column 0 -pady 2 -padx 2
  grid $helpers.reps -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid $helpers.refresh_reps -row $gridrow -column 2 -pady 2 -padx 2
  incr gridrow

  # Populate initial list of selection texts from reps
  refresh_reps

  ############# Atoms from atom, bond, angle, dihedral labels ####################
  ttk::button $helpers.labeled_atoms -text "Insert labeled atoms" -command {::cv_dashboard::insert_labels Atoms}
  ttk::button $helpers.labeled_var -text "Insert labeled..." -command {::cv_dashboard::insert_labels combo}
  ttk::combobox $helpers.labels -justify left -state readonly
  $helpers.labels configure -values [list Bonds Angles Dihedrals]
  $helpers.labels set Bonds

  grid $helpers.labeled_atoms -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid $helpers.labeled_var -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid $helpers.labels -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew
  incr gridrow

  ################# Insert file name from file picker ###########################
  ttk::radiobutton $helpers.files1 -variable ::cv_dashboard::filetype -text "atomsFile" -value "atomsFile"
  ttk::radiobutton $helpers.files2 -variable ::cv_dashboard::filetype -text "refPositionsFile" -value "refPositionsFile"
  ttk::button $helpers.insert_file -text "Pick file" \
    -command [list ::cv_dashboard::insert_filename] -padding "2 0 2 0"

  grid $helpers.files1 -row $gridrow -column 0 -pady 2 -padx 2
  grid $helpers.files2 -row $gridrow -column 1 -pady 2 -padx 2
  grid $helpers.insert_file -row $gridrow -column 2 -pady 2 -padx 2
  incr gridrow

  grid columnconfigure $helpers 0 -weight 1
  grid columnconfigure $helpers 1 -weight 1
  grid columnconfigure $helpers 2 -weight 1

  grid $docs -sticky ew
  grid $helpers -sticky ew
  grid columnconfigure $w.editor.fl 0 -weight 1


  # Right frame: text widget w scrollbar and Apply/Cancel buttons
  frame $w.editor.fr
  tk::text $w.editor.fr.text -undo 1 -yscrollcommand [list $w.editor.fr.vsb set] -background white
  ttk::scrollbar $w.editor.fr.vsb -orient vertical -command [list $w.editor.fr.text yview]
  $w.editor.fr.text insert 1.0 $cfg
  set ::cv_dashboard::being_edited $cvs
  grid $w.editor.fr.text -row 0 -columnspan 3 -sticky nsew
  grid $w.editor.fr.vsb -row 0 -column 3 -sticky nsew

  # Ctrl-s anywhere in the window saves/applies
  bind $w.editor <Control-s> ::cv_dashboard::edit_apply
  bind $w.editor <Control-a> "$w.editor.fr.text tag add sel 1.0 end-1c"

  set gridrow 1
  ttk::button $w.editor.fr.apply -text "Apply \[Ctrl-s\]" -command ::cv_dashboard::edit_apply -padding "2 0 2 0"
  ttk::button $w.editor.fr.cancel -text "Cancel" -command ::cv_dashboard::edit_cancel -padding "2 0 2 0"
  ttk::button $w.editor.fr.clear -text "Clear" -command "$w.editor.fr.text delete 1.0 end" -padding "2 0 2 0"

  grid $w.editor.fr.apply -row $gridrow -column 0 -sticky e -pady 2 -padx 2
  grid $w.editor.fr.cancel -row $gridrow -column 1 -sticky w -pady 2 -padx 2
  grid $w.editor.fr.clear -row $gridrow -column 2 -sticky w -pady 2 -padx 2

  grid columnconfigure $w.editor.fr 0 -weight 1
  grid columnconfigure $w.editor.fr 1 -weight 1
  grid columnconfigure $w.editor.fr 2 -weight 1
  grid rowconfigure $w.editor.fr 0 -weight 1

  pack $w.editor.fl -fill both -side left
  pack $w.editor.fr -fill both -expand yes -padx 2 -pady 2
  # pack $w.editor.fr -side bottom -fill both -expand yes
}


# Multi-platform solution from http://wiki.tcl.tk/557
proc ::cv_dashboard::invokeBrowser {url} {
  # open is the OS X equivalent to xdg-open on Linux, start is used on Windows
  set commands {xdg-open open start}
  foreach browser $commands {
    if {$browser eq "start"} {
      set command [list {*}[auto_execok start] {}]
    } else {
      set command [auto_execok $browser]
    }
    if {[string length $command]} {
      break
    }
  }

  if {[string length $command] == 0} {
    return -code error "couldn't find browser"
  }
  if {[catch {exec {*}$command $url &} error]} {
    return -code error "couldn't execute '$command': $error"
  }
}


# Insert atomNumbers command for given selection text
proc ::cv_dashboard::atoms_from_sel { source } {
  set w .cv_dashboard_window

  # Called from textbox
  if { $source == "textbox" } {
    set seltext [$w.editor.fl.helpers.seltext get]
  } elseif { $source == "reps" } {
    set seltext [$w.editor.fl.helpers.reps get]
  }

  if {[llength $seltext] == 0 } {
    return
  }
  set sel [atomselect top $seltext]
  set serials [$sel get serial]
  $sel delete

  if {[llength $serials] == 0 } {
    tk_messageBox -icon error -title "Colvars error" -parent .cv_dashboard_window\
      -message "Selection text matches zero atoms"
    return
  }
  editor_replace "        # Selection: \"$seltext\"\n        atomNumbers $serials\n"
}


# Replace selection in editor with given string
proc ::cv_dashboard::editor_replace { text } {
  set t .cv_dashboard_window.editor.fr.text

  if {[$t tag ranges sel] != ""} {
    $t delete sel.first sel.last
  }
  $t insert insert $text
}


# Insert atom numbers or components from currently labeled objects in VMD
proc ::cv_dashboard::insert_labels {obj} {
  set w .cv_dashboard_window
  if {$obj == "combo"} {
    set obj [$w.editor.fl.helpers.labels get]
  }

  if { $obj == "Atoms" } {
    set serials [list]
    foreach l [label list $obj] {
      if { [lindex $l 2] == "hide" } { continue }
      set a [lindex $l 0]
      lappend serials [expr [lindex $a 1] + 1] ;# going from VMD 0-based to 1-based atomNumbers
    }
    if {[llength $serials] > 0} {
      editor_replace "        # Atom labels\n        atomNumbers $serials\n"
    }
  } else {
    set n(Bonds) 2
    set n(Angles) 3
    set n(Dihedrals) 4
    set cvc(Bonds) distance
    set cvc(Angles) angle
    set cvc(Dihedrals) dihedral
    foreach l [label list $obj] {
      if { [lindex $l 2] == "hide" } { continue }
      set cfg "    $cvc($obj) \{\n"
      for {set i 0} { $i < $n($obj) } {incr i} {
        set a [lindex $l $i]
        set serial [expr [lindex $a 1] + 1]
        append cfg "        group[expr $i+1] \{\n            atomNumbers $serial\n        \}\n"
      }
      append cfg "    \}\n"
      editor_replace $cfg
    }
  }
}

# Insert contents of template file
proc ::cv_dashboard::insert_template {} {
  set w .cv_dashboard_window
  if { [info exists ::cv_dashboard::template_dir] } {
    set path [tk_getOpenFile -initialdir $::cv_dashboard::template_dir]
  } else {
    set path [tk_getOpenFile -initialdir [pwd]]
  }
  if [string compare $path ""] {
    # Save directory for next invocation of this dialog
    set ::cv_dashboard::template_dir [file dirname $path]
    set in [open $path r]
    editor_replace [read $in]
    close $in
  }
}


# Insert filename
proc ::cv_dashboard::insert_filename {} {
  variable ::cv_dashboard::filetype
  set w .cv_dashboard_window

  if { [info exists ::cv_dashboard::atomfile_dir] } {
    set path [tk_getOpenFile -filetypes {{"PDB" .pdb} {"XYZ" .xyz} {"All files" *}} \
        -initialdir $::cv_dashboard::atomfile_dir]
  } else {
    set path [tk_getOpenFile -filetypes {{"PDB" .pdb} {"XYZ" .xyz} {"All files" *}} -initialdir [pwd]]
  }
  if [string compare $path ""] {
    # Save directory for next invocation of this dialog
    set ::cv_dashboard::atomfile_dir [file dirname $path]
    set coltype [string range $filetype 0 end-4]
    editor_replace "    $filetype $path\n    ${coltype}Col O\n    ${coltype}ColValue 1\n"
  }
}


# Submit current config text from editor to cvm
proc ::cv_dashboard::edit_apply {} {
  set w .cv_dashboard_window
  foreach c $::cv_dashboard::being_edited {
    run_cv colvar $c delete
  }
  set cfg [$w.editor.fr.text get 1.0 end-1c]
  if { $cfg != "" } {
    set res [run_cv config $cfg]
    if { [string compare $res ""] } {
      # error: restore backed up cfg
      run_cv config $::cv_dashboard::backup_cfg
      refresh_table
      # Do not destroy editor window (give user a chance to fix input)
      return
    }
  }
  set ::cv_dashboard::being_edited {}
  destroy $w.editor
  refresh_table
}


# Close editor without applying
proc ::cv_dashboard::edit_cancel {} {
  set w .cv_dashboard_window
  set ::cv_dashboard::being_edited {}
  destroy $w.editor
}


proc ::cv_dashboard::refresh_reps {} {
  set w .cv_dashboard_window
  set numreps [molinfo top get numreps]
  set reps [list]
  for {set i 0} {$i < $numreps} {incr i} {
    lappend reps [lindex [molinfo top get [list [list selection $i]]] 0]
  }
  $w.editor.fl.helpers.reps configure -values $reps
}

#################################################################
# Interactive plot window
#################################################################


# Create plot window
proc ::cv_dashboard::plot { { type timeline } } {
  variable ::cv_dashboard::plothandle
  set ::cv_dashboard::plottype $type

  # Remove existing plot, if any
  if { [info exists plothandle] } {
    catch {$plothandle quit}
    unset plothandle
  }

  set cvs [selected_colvars]
  if { [llength $cvs] == 0 } {
    # If no selection, plot all variables
    set cvs $::cv_dashboard::cvs
    if { [llength $cvs] == 0 } {
      return
    }
  }

  # Analyze colvar values to split vector values into scalars with numeric index
  # store array of names for each scalar value
  set total_dim 0
  set name_list {}
  foreach c $cvs {
    set val [run_cv colvar $c update]
    set comps($c) [selected_comps $c]
    set size [llength $comps($c)]
    incr total_dim $size
    if { $size == 1 } {
      set names($c) $c
      lappend name_list $c
      set y($c) {}
    } else {
      set names($c) {}
      for {set i 0} {$i < $size} {incr i} {
        set n "${c}_[expr [lindex $comps($c) $i] + 1]"
        lappend names($c) $n
        lappend name_list $n
        set y($n) {}
      }
    }
  }

  if { $type == "2cv" && $total_dim != 2 } {
    tk_messageBox -icon error -title "Colvars Dashboard Error"\
      -message "Select exactly 2 scalar quantities for pairwise plot.\n"
    return
  }

  set nf [molinfo top get numframes]
  # Get list of values for all frames
  for {set f 0} {$f< $nf} {incr f} {
    run_cv frame $f
    foreach c $cvs {
      set val [run_cv colvar $c update]
      foreach ni $names($c) vi $comps($c) {
        lappend y($ni) [lindex $val $vi]
      }
    }
  }

  if { $type == "timeline"} {
    set plothandle [multiplot \
      -title {Colvars trajectory   [click, keyb arrows (+ Shift/Ctrl) to navigate & zoom, v/h to fit vert/horizontally]} \
      -xlabel "Frame" -ylabel "Value" -nostats]
    set x {}
    for {set f 0} {$f < $nf} {incr f} { lappend x $f }
    foreach n $name_list {
      $plothandle add $x $y($n) -legend $n
    }
  } elseif { $type == "2cv"} {
    set xname [lindex $name_list 0]
    set yname [lindex $name_list 1]
    set plothandle [multiplot -title {Colvars trajectory   [click on markers, keyb arrows (+ Shift/Ctrl) to navigate]} \
    -xlabel $xname -ylabel $yname -nostats -marker circle -fill white -radius 4 -callback ::cv_dashboard::marker_clicked]
    $plothandle add $y($xname) $y($yname)
  }

  $plothandle replot
  # bind mouse and keyboard events to callbacks
  set plot_ns [namespace qualifiers $::cv_dashboard::plothandle]

  if { $type == "timeline" } {
    bind [set ${plot_ns}::w] <Button-1>       { ::cv_dashboard::plot_clicked %x %y }
    bind [set ${plot_ns}::w] <Left>           { ::cv_dashboard::chg_frame -1 }
    bind [set ${plot_ns}::w] <Right>          { ::cv_dashboard::chg_frame 1 }
    bind [set ${plot_ns}::w] <Shift-Left>     { ::cv_dashboard::chg_frame -10 }
    bind [set ${plot_ns}::w] <Shift-Right>    { ::cv_dashboard::chg_frame 10 }
    bind [set ${plot_ns}::w] <Control-Left>   { ::cv_dashboard::chg_frame -50 }
    bind [set ${plot_ns}::w] <Control-Right>  { ::cv_dashboard::chg_frame 50 }
    bind [set ${plot_ns}::w] <Home>           { ::cv_dashboard::chg_frame start }
    bind [set ${plot_ns}::w] <End>            { ::cv_dashboard::chg_frame end }
    bind [set ${plot_ns}::w] <Up>             { ::cv_dashboard::zoom 0.25 }
    bind [set ${plot_ns}::w] <Down>           { ::cv_dashboard::zoom 4 }
    bind [set ${plot_ns}::w] <Shift-Up>       { ::cv_dashboard::zoom 0.0625 }
    bind [set ${plot_ns}::w] <Shift-Down>     { ::cv_dashboard::zoom 16 }
    bind [set ${plot_ns}::w] <v>              { ::cv_dashboard::fit_vertically }
    bind [set ${plot_ns}::w] <h>              { ::cv_dashboard::fit_horizontally }
  } elseif { $type == "2cv"} {
    bind [set ${plot_ns}::w] <Left>           { ::cv_dashboard::chg_frame -1 }
    bind [set ${plot_ns}::w] <Right>          { ::cv_dashboard::chg_frame 1 }
    bind [set ${plot_ns}::w] <Shift-Left>     { ::cv_dashboard::chg_frame -10 }
    bind [set ${plot_ns}::w] <Shift-Right>    { ::cv_dashboard::chg_frame 10 }
    bind [set ${plot_ns}::w] <Control-Left>   { ::cv_dashboard::chg_frame -50 }
    bind [set ${plot_ns}::w] <Control-Right>  { ::cv_dashboard::chg_frame 50 }
    bind [set ${plot_ns}::w] <Home>           { ::cv_dashboard::chg_frame start }
    bind [set ${plot_ns}::w] <End>            { ::cv_dashboard::chg_frame end }
  }

  # Update frame to display frame marker in new plot
  update_frame internal [molinfo top] w
}


# Callback for click inside plot window, at coords x y
proc ::cv_dashboard::plot_clicked { x y } {

  set ns [namespace qualifiers $::cv_dashboard::plothandle]
  set xplotmin [set ${ns}::xplotmin]
  set xplotmax [set ${ns}::xplotmax]
  set yplotmin [set ${ns}::yplotmin]
  set yplotmax [set ${ns}::yplotmax]
  set scalex [set ${ns}::scalex]
  set xmin [set ${ns}::xmin]

  # note: ymax < ymin because of pixel coordinate convention
  if { [expr {($x < $xplotmin) || ($x > $xplotmax) || ($y > $yplotmin) || ($y < $yplotmax)}] } {
    return
  }

  animate goto [expr { ($x - $xplotmin) / $scalex + $xmin}]
  if { $::cv_dashboard::track_frame == 0 } {
    # frame change doesn't trigger refresh, so we refresh manually
    refresh_table
  }
}


# Callback for click on marker in 2 cv plot
proc ::cv_dashboard::marker_clicked { index x y color marker } {

  animate goto [expr {$index - 1 }]
  if { $::cv_dashboard::track_frame == 0 } {
    # frame change doesn't trigger refresh, so we refresh manually
    refresh_table
  }
}


# Change frame in reaction to user input (arrow keys)
proc ::cv_dashboard::chg_frame { shift } {

  set nf [molinfo top get numframes]

  if { $shift == "start" } {
    set f 0
  } elseif { $shift == "end" } {
    set f [expr $nf - 1]
  } else {
    set f [expr $::cv_dashboard::current_frame + $shift]
  }

  # Keep within bounds [[O, numframes[[
  if { $f < 0 } { set f 0 }
  if { $f >= $nf } { set f [expr $nf - 1] }

  animate goto $f
  if { $::cv_dashboard::track_frame == 0 } {
    # frame change doesn't trigger refresh, so we refresh manually
    refresh_table
  }
}


# Change zoom in reaction to user input; recenter marker
proc ::cv_dashboard::zoom { factor } {
  variable ::cv_dashboard::plothandle

  set ns [namespace qualifiers $plothandle]
  set xmin [set ${ns}::xmin]
  set xmax [set ${ns}::xmax]

  set f $::cv_dashboard::current_frame
  # rescale current half-width
  set half_width [expr { int (($xmax - $xmin) * $factor / 2)}]
  if { $half_width < 1 } {
    # Don't collapse the x axis to a single point
    return
  }
  set fmin [expr { $f - $half_width }]
  if {$fmin < 0} { set fmin 0 }

  set fmax [expr { $f + $half_width }]
  set max_f [expr [molinfo top get numframes] - 1]
  if {$fmax > $max_f} { set fmax $max_f }

  $plothandle configure -xmin $fmin -xmax $fmax -plot
  update_frame internal [molinfo top] w
}


# Fit vertical axis to y values within the current x range
proc ::cv_dashboard::fit_vertically {} {
  variable ::cv_dashboard::plothandle

  set ns [namespace qualifiers $plothandle]
  set xmin [set ${ns}::xmin]
  set xmax [set ${ns}::xmax]
  set ydata [$plothandle ydata]
  set ymin [lindex [lindex $ydata 0] $xmin]
  set ymax $ymin
  foreach yset $ydata {
    for { set f $xmin } { $f <= $xmax } { incr f } {
      set y [lindex $yset $f]
      if { $y > $ymax } { set ymax $y }
      if { $y < $ymin } { set ymin $y }
    }
  }
  $plothandle configure -ymin $ymin -ymax $ymax -plot
  update_frame internal [molinfo top] w
}


# Fit all x values within horizontal range, then fit vertically
proc ::cv_dashboard::fit_horizontally {} {
  variable ::cv_dashboard::plothandle

  $plothandle configure -xmin auto -xmax auto -ymin auto -ymax auto -plot
  update_frame internal [molinfo top] w
}


# Display frame marker in plot at given frame
proc ::cv_dashboard::display_marker { f } {
  variable ::cv_dashboard::plothandle
  if [info exists plothandle] {
    # detect if plot was closed
    if [catch {$plothandle getpath}] {
      unset plothandle
    } else {
      # we tinker a little with Multiplot's internals to get access to its Tk canvas
      # necessary because Multiplot does not expose an interface to draw & delete
      # objects without redrawing the whole plot - which takes too long for this
      set ns [namespace qualifiers $plothandle]
      if { $::cv_dashboard::plottype == "timeline" } {

        set xmin [set ${ns}::xmin]
        set xmax [set ${ns}::xmax]
        # Move plot boundaries if necessary
        if { $f < $xmin } {
          set xmax [expr { $xmax + $f - $xmin }]
          set xmin $f
          $plothandle configure -xmin $xmin -xmax $xmax -plot
        }
        if { $f > $xmax } {
          set xmin [expr { $xmin + $f - $xmax }]
          set xmax $f
          $plothandle configure -xmin $xmin -xmax $xmax -plot
        }

        set y1 [set ${ns}::yplotmin]
        set y2 [set ${ns}::yplotmax]
        set xplotmin [set ${ns}::xplotmin]
        set scalex [set ${ns}::scalex]
        set x [expr $xplotmin+($scalex*($f-$xmin))]

        set canv "[set ${ns}::w].f.cf"
        $canv delete frame_marker
        $canv create line  $x $y1 $x $y2 -fill blue -tags frame_marker
      } elseif { $::cv_dashboard::plottype == "2cv" } {
        set x [lindex [ lindex [$plothandle xdata] 0] $f]
        set y [lindex [ lindex [$plothandle ydata] 0] $f]

        set xmin [set ${ns}::xmin]
        set ymin [set ${ns}::ymin]

        set radius 5
        set xplotmin [set ${ns}::xplotmin]
        set scalex [set ${ns}::scalex]
        set x1 [expr {$xplotmin+$scalex*($x-$xmin) - $radius}]
        set x2 [expr {$xplotmin+$scalex*($x-$xmin) + $radius}]

        set yplotmin [set ${ns}::yplotmin]
        set scaley [set ${ns}::scaley]
        set y1 [expr {$yplotmin+$scaley*($y-$ymin) - $radius}]
        set y2 [expr {$yplotmin+$scaley*($y-$ymin) + $radius}]

        set canv "[set ${ns}::w].f.cf"
        $canv delete frame_marker
        $canv create oval $x1 $y1 $x2 $y2 -outline white -fill blue -tags frame_marker
      }
    }
  }
}

