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

# TODO:
# - histograms
# - graphical representations such as rotation_display
# - show atom groups as representations

package provide cv_dashboard 1.0


namespace eval ::cv_dashboard {
  # General UI state
  variable current_frame 0  ;# linked to frame display
  variable cv_values {}     ;# linked to value table
  variable cvs {}           ;# linked to colvar table
  variable track_frame 1    ;# start tracking by default

  # State variables for config editor
  variable being_edited
  variable backup_cfg
  variable filetype "atomsFile"

  # Handle to keep track of a single interactive plot
  variable plothandle
  variable plottype         ;# either timeline or 2cv
}


#################################################################
# Main UI: the dashboard
#################################################################


proc cv_dashboard {} {
  # If window already exists, destroy it
  catch { destroy .cv_dashboard_window }
  ::cv_dashboard::createWindow
}


# Creat main window
proc ::cv_dashboard::createWindow {} {

  package require tablelist

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
  grid [ttk::button $w.load -text "Load config file" -command ::cv_dashboard::load -padding "2 0 2 0"] \
    -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $w.save -text "Save colvars config" -command ::cv_dashboard::save -padding "2 0 2 0"] \
    -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $w.reset -text "Reset Colvars Module" -command ::cv_dashboard::reset -padding "2 0 2 0"] \
    -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew

  tablelist::tablelist $w.cvtable -columns {
    0 Colvars
  } -stretch all -selectmode extended -selecttype row -listvariable ::cv_dashboard::cvs
  tablelist::tablelist $w.valuetable -columns {
    0 Values
  } -stretch all -listvariable ::cv_dashboard::cv_values

  refresh_table
  bind [$w.cvtable bodytag] <Double-1> ::cv_dashboard::edit

  incr gridrow
  grid $w.cvtable -row $gridrow -column 0 -sticky news
  grid $w.valuetable -row $gridrow -column 1 -sticky news -columnspan 2
  grid rowconfigure $w $gridrow -weight 1

  incr gridrow
  grid [ttk::button $w.edit -text "Edit \[dbl-click\]" -command ::cv_dashboard::edit -padding "2 0 2 0"] \
    -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $w.add -text "New" -command ::cv_dashboard::add -padding "2 0 2 0"] \
    -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $w.del -text "Delete" -command ::cv_dashboard::del -padding "2 0 2 0"] \
    -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew

  incr gridrow
  grid [ttk::button $w.plot -text "Timeline plot" -command ::cv_dashboard::plot -padding "2 0 2 0"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $w.plot2cv -text "2d plot" -command {::cv_dashboard::plot 2cv} -padding "2 0 2 0"] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [ttk::button $w.refresh -text "Refresh table" -command ::cv_dashboard::refresh_table -padding "2 0 2 0"] -row $gridrow -column 2 -pady 2 -padx 2 -sticky nsew

  incr gridrow
  grid [label $w.frameTxt -text "Frame:"] -row $gridrow -column 0 -pady 2 -padx 2 -sticky nsew
  grid [label $w.frame -textvariable ::cv_dashboard::current_frame] -row $gridrow -column 1 -pady 2 -padx 2 -sticky nsew
  grid [ttk::checkbutton $w.trackFrame -text "Track" -command ::cv_dashboard::change_track_frame -variable ::cv_dashboard::track_frame] \
    -row $gridrow -column 2  -pady 2 -padx 2 -sticky nsew
  change_track_frame ;# activate tracking if necessary

  grid columnconfigure $w 0 -weight 1
  grid columnconfigure $w 1 -weight 1
  grid columnconfigure $w 2 -weight 1
}


# Refresh the table with a list of existing CVs and their values
proc ::cv_dashboard::refresh_table {} {
  set w .cv_dashboard_window
  set ::cv_dashboard::cvs [run_cv list]
  update_frame vmd_frame [molinfo top] w
}


# Refresh the table with a list of colvar values
proc ::cv_dashboard::refresh_values {} {
  run_cv update
  set w .cv_dashboard_window
  set ::cv_dashboard::cv_values {}
  foreach c [run_cv list] {
    lappend ::cv_dashboard::cv_values [format_value [run_cv colvar $c value] ]
  }
}


# Format colvar values returned by cv for display in table
proc ::cv_dashboard::format_value val {
  if {[llength $val] == 1} {
    return [format "%.4g" $val]
  } else {
    set s "{("
    foreach v [lrange $val 0 end-1] {
      append s "[format_value $v] "
    }
    append s "[format_value [lindex $val end]])}"
    return $s
  }
}


# Return list of selected colvars in the table
proc ::cv_dashboard::selected {} {
  set w .cv_dashboard_window
  set l {}
  refresh_table ;# to make sure selected colvars haven't been deleted
  foreach i [$w.cvtable curselection] {
    lappend l [lindex $::cv_dashboard::cvs $i]
  }
  return $l
}


# Enable or disable real-time tracking of VMD frame
proc ::cv_dashboard::change_track_frame {} {
  global vmd_frame

  set molid [molinfo top]
  if {$::cv_dashboard::track_frame} {
    trace add variable vmd_frame($molid) write ::cv_dashboard::update_frame
    update_frame vmd_frame [molinfo top] w
  } else {
    trace remove variable vmd_frame($molid) write ::cv_dashboard::update_frame
  }
}


# Load config from file
proc ::cv_dashboard::load {} {
  set path [tk_getOpenFile -filetypes {{"Colvars cfg" .in} {"Colvars cfg" .colvars} {"All files" *}}]
  if [string compare $path ""] {
    run_cv configfile $path
    refresh_table
  }
}


# Save config of colvars to file (can we do it w/ biases ? need bias type keyword)
proc ::cv_dashboard::save {} {

  set path [tk_getSaveFile -filetypes {{"Colvars cfg" .in} {"Colvars cfg" .colvars} {"All files" *}}]

  if [string compare $path ""] {
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
  foreach c [selected] {
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


# Reset cvm
proc ::cv_dashboard::reset {} {
  run_cv reset
  refresh_table
}


#################################################################
# General utilities
#################################################################


# Call the "cv" interface to Colvars, catching errors and displaying them to the user
proc run_cv args  {
  if [ catch { cv {*}$args } res ] {
    tk_messageBox -icon error -title "Colvars error" -parent .cv_dashboard_window\
      -message "Error running command:\n\n$args" -detail "$res"
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
    # Provide simple template
    set cfg "# You can edit or replace the example colvar config below.\n\
colvar {\n  name d\n  distance {\n    group1 { atomNumbers 1 2 }\n    group2 { atomNumbers 3 4 }\n  }\n}\n"
  } else {
    set cvs [selected]
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
  set editor [toplevel $w.editor]
  wm title $editor "Colvar config editor"

  # Left frame: utility buttons
  frame $w.editor.fl
  set gridrow 0

  ttk::button $w.editor.fl.onlinedoc1 -text "Online doc: defining collective variables" -padding "4 2 4 2"\
    -command [list ::cv_dashboard::invokeBrowser "http://colvars.github.io/colvars-refman-vmd/colvars-refman-vmd.html#sec:colvar"]
  ttk::button $w.editor.fl.onlinedoc2 -text "Online doc: defining atom groups" -padding "4 2 4 2"\
    -command [list ::cv_dashboard::invokeBrowser "http://colvars.github.io/colvars-refman-vmd/colvars-refman-vmd.html#sec:colvar_atom_groups"]
  ttk::button $w.editor.fl.onlinedoc3 -text "Online doc: types of variables (components)" -padding "4 2 4 2"\
    -command [list ::cv_dashboard::invokeBrowser "http://colvars.github.io/colvars-refman-vmd/colvars-refman-vmd.html#sec:cvc"]

  grid $w.editor.fl.onlinedoc1 -row $gridrow -column 0 -columnspan 3 -pady 5
  incr gridrow
  grid $w.editor.fl.onlinedoc2 -row $gridrow -column 0 -columnspan 3 -pady 5
  incr gridrow
  grid $w.editor.fl.onlinedoc3 -row $gridrow -column 0 -columnspan 3 -pady 5
  incr gridrow


  tk::label $w.editor.fl.seltext_label -text "Selection text:"
  tk::entry $w.editor.fl.seltext -bg white
  # Bind Return key in seltext entry to proc creating the atomNumbers line
  bind $w.editor.fl.seltext <Return> "::cv_dashboard::atoms_from_sel textbox"
  ttk::button $w.editor.fl.fromsel -text "Insert atoms \[Enter\]" \
    -command "::cv_dashboard::atoms_from_sel textbox" -padding "2 0 2 0"

  grid $w.editor.fl.seltext_label -row $gridrow -column 0 -pady 2 -padx 2
  grid $w.editor.fl.seltext -row $gridrow -column 1 -sticky ew -pady 2 -padx 2
  grid $w.editor.fl.fromsel -row $gridrow -column 2 -pady 2 -padx 2
  incr gridrow


  # Double click or Enter to insert
  ttk::button $w.editor.fl.refresh_reps -text "Refresh list" -command ::cv_dashboard::refresh_reps
  tk::listbox $w.editor.fl.reps -selectmode browse
  bind $w.editor.fl.reps <Return> "::cv_dashboard::atoms_from_sel reps"
  bind $w.editor.fl.reps <Double-1> "::cv_dashboard::atoms_from_sel reps"

  grid $w.editor.fl.refresh_reps -row $gridrow -column 0 -pady 2 -padx 2
  grid $w.editor.fl.reps -row $gridrow -column 1 -columnspan 2 -pady 2 -padx 2 -sticky nsew
  incr gridrow

  # Populate initial list of selection texts from reps
  refresh_reps

  ttk::radiobutton $w.editor.fl.files1 -variable ::cv_dashboard::filetype -text "atomsFile" -value "atomsFile"
  ttk::radiobutton $w.editor.fl.files2 -variable ::cv_dashboard::filetype -text "refPositionsFile" -value "refPositionsFile"
  ttk::button $w.editor.fl.insert_file -text "Pick file" \
    -command [list ::cv_dashboard::insert_filename] -padding "2 0 2 0"

  grid $w.editor.fl.files1 -row $gridrow -column 0 -pady 2 -padx 2
  grid $w.editor.fl.files2 -row $gridrow -column 1 -pady 2 -padx 2
  grid $w.editor.fl.insert_file -row $gridrow -column 2 -pady 2 -padx 2
  incr gridrow


  grid columnconfigure $w.editor.fl 0 -weight 1
  grid columnconfigure $w.editor.fl 1 -weight 1
  grid columnconfigure $w.editor.fl 2 -weight 1


  # Right frame: text widget w scrollbar and Apply/Cancel buttons
  frame $w.editor.fr
  tk::text $w.editor.fr.text -undo 1 -yscrollcommand [list $w.editor.fr.vsb set] -background white
  ttk::scrollbar $w.editor.fr.vsb -orient vertical -command [list $w.editor.fr.text yview]
  $w.editor.fr.text insert 1.0 $cfg
  set ::cv_dashboard::being_edited $cvs
  grid $w.editor.fr.text -row 0 -columnspan 2 -sticky nsew
  grid $w.editor.fr.vsb -row 0 -column 2 -sticky nsew

  # Ctrl-s anywhere in the window saves/applies
  bind $w.editor <Control-s> ::cv_dashboard::edit_apply

  set gridrow 1
  ttk::button $w.editor.fr.apply -text "Apply \[Ctrl-s\]" -command ::cv_dashboard::edit_apply -padding "2 0 2 0"
  ttk::button $w.editor.fr.cancel -text "Cancel" -command ::cv_dashboard::edit_cancel -padding "2 0 2 0"
  grid $w.editor.fr.apply -row $gridrow -column 0 -sticky e -pady 2 -padx 2
  grid $w.editor.fr.cancel -row $gridrow -column 1 -sticky w -pady 2 -padx 2

  grid columnconfigure $w.editor.fr 0 -weight 1
  grid columnconfigure $w.editor.fr 1 -weight 1
  grid rowconfigure $w.editor.fr 0 -weight 1

  pack $w.editor.fl -fill both -side left
  pack $w.editor.fr -fill both -side left -expand yes
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
  puts "$command $url"
  if {[catch {exec {*}$command $url &} error]} {
    return -code error "couldn't execute '$command': $error"
  }
}


# Insert atomNumbers command for given selection text
proc ::cv_dashboard::atoms_from_sel { source } {
  set w .cv_dashboard_window

  # Called from textbox
  if { $source == "textbox" } {
    set seltext [$w.editor.fl.seltext get]
  } elseif { $source == "reps" } {
    set seltext [$w.editor.fl.reps get active]
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
  $w.editor.fr.text insert insert "      # $seltext\n      atomNumbers $serials\n"
}


# Insert filename
proc ::cv_dashboard::insert_filename {} {
  variable ::cv_dashboard::filetype
  set w .cv_dashboard_window
  set path [tk_getOpenFile -filetypes {{"PDB" .pdb} {"All files" *}}]
  if { ![string compare $path ""] } {
    return
  }
  set coltype [string range $filetype 0 end-4]
  $w.editor.fr.text insert insert "    $filetype $path\n    ${coltype}Col O\n    ${coltype}ColValue 1\n"
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
  $w.editor.fl.reps delete 0 end
  $w.editor.fl.reps insert 0 {*}$reps
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

  set cvs [selected]
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
    set size [llength $val]
    incr total_dim $size
    if { $size == 1 } {
      set names($c) $c
      lappend name_list $c
      set y($c) {}
    } else {
      for {set i 1} {$i <= $size} {incr i} {
        set n $c; append n "_$i"
        lappend names($c) $n
        lappend name_list $n
        set y($n) {}
      }
    }
  }

  if { $type == "2cv" } {
    if { $total_dim > 2 } {
      puts "Warning: 2d plot will use the first 2 of $total_dim scalar dimensions in selection."
    } elseif { $total_dim < 2 } {
      tk_messageBox -icon error -title "Colvars Dashboard Error"\
        -message "Not enough data selected. 2 scalar sets are needed for 2d plot.\n"
      return
    }
  }

  set nf [molinfo top get numframes]
  # Get list of values for all frames
  for {set f 0} {$f< $nf} {incr f} {
    run_cv frame $f
    foreach c $cvs {
      set val [run_cv colvar $c update]
      foreach ni $names($c) vi $val {
        lappend y($ni) $vi
      }
    }
  }

  if { $type == "timeline"} {
    set plothandle [multiplot \
      -title {Colvars trajectory   [left-click, keyb arrows (+ Shift/Ctrl) to navigate & zoom, v/h to fit vert/horizontally]} \
      -xlabel "Frame" -ylabel "Value" -nostats]
    set x {}
    for {set f 0} {$f < $nf} {incr f} { lappend x $f }
    foreach n $name_list {
      $plothandle add $x $y($n) -legend $n
    }
  } elseif { $type == "2cv"} {
    set xname [lindex $name_list 0]
    set yname [lindex $name_list 1]
    set plothandle [multiplot -title {Colvars trajectory} \
    -xlabel $xname -ylabel $yname -nostats -marker circle]
    $plothandle add $y($xname) $y($yname)
  }

  $plothandle replot

  if { $type == "timeline" } {
    # bind mouse and keyboard events to callbacks
    set plot_ns [namespace qualifiers $::cv_dashboard::plothandle]
    bind [set ${plot_ns}::w] <Button-1>       { ::cv_dashboard::plot_clicked %x %y }
    bind [set ${plot_ns}::w] <Left>           { ::cv_dashboard::chg_frame -1 }
    bind [set ${plot_ns}::w] <Right>          { ::cv_dashboard::chg_frame 1 }
    bind [set ${plot_ns}::w] <Shift-Left>     { ::cv_dashboard::chg_frame -10 }
    bind [set ${plot_ns}::w] <Shift-Right>    { ::cv_dashboard::chg_frame 10 }
    bind [set ${plot_ns}::w] <Control-Left>   { ::cv_dashboard::chg_frame -50 }
    bind [set ${plot_ns}::w] <Control-Right>  { ::cv_dashboard::chg_frame 50 }
    bind [set ${plot_ns}::w] <Up>             { ::cv_dashboard::zoom 0.25 }
    bind [set ${plot_ns}::w] <Down>           { ::cv_dashboard::zoom 4 }
    bind [set ${plot_ns}::w] <Shift-Up>       { ::cv_dashboard::zoom 0.0625 }
    bind [set ${plot_ns}::w] <Shift-Down>     { ::cv_dashboard::zoom 16 }
    bind [set ${plot_ns}::w] <v>              { ::cv_dashboard::fit_vertically }
    bind [set ${plot_ns}::w] <h>              { ::cv_dashboard::fit_horizontally }
  }
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


# Change frame in reaction to user input (arrow keys)
proc ::cv_dashboard::chg_frame { shift } {
  set f [expr $::cv_dashboard::current_frame + $shift]

  # Keep within bounds [[O, numframes[[
  if { $f < 0 } { set f 0 }
  set nf [molinfo top get numframes]
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
  update_frame vmd_frame [molinfo top] w
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
  update_frame vmd_frame [molinfo top] w
}


# Fit all x values within horizontal range, then fit vertically
proc ::cv_dashboard::fit_horizontally {} {
  variable ::cv_dashboard::plothandle

  $plothandle configure -xmin auto -xmax auto -ymin auto -ymax auto -plot
  update_frame vmd_frame [molinfo top] w
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
#         set xmin [set ${ns}::xmin]
#         set xmax [set ${ns}::xmax]
#         set ymin [set ${ns}::ymin]
#         set ymax [set ${ns}::ymax]
#         set x [lindex set ${ns}::]
#         set y2 [set ${ns}::yplotmax]
      }
    }
  }
} 

# Create a window to begin with - until we have a menu entry
cv_dashboard
