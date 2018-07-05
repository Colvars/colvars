# GUI for the Colvars Module in VMD
# Jérôme Hénin <henin@ibpc.fr> 2018

# At this stage, only acts on the "top" molecule
# could query cvm to know what molecule is concerned

namespace eval ::cvgui {
  variable current_frame 0
  variable cv_values {}
  variable cvs {}
  variable track_frame 1 ;# start tracking by default
  variable being_edited
  variable backup_cfg
  variable plothandle
}


# Call the "cv" interface to Colvars, catching errors and displaying them to the user
proc run_cv args  {
  if [ catch { cv {*}$args } res ] {
    tk_messageBox -icon error -title "Error" -parent .cvgui_window\
      -message "Error running command:\n\n$args\n\n$res"
    return -1
  }
  return $res
}


# Refresh the table with a list of existing CVs and their values
proc ::cvgui::refresh_table {} {
  set w .cvgui_window
  set ::cvgui::cvs [run_cv list]
  update_frame vmd_frame [molinfo top] w
}

proc ::cvgui::refresh_values {} {
  run_cv update
  set w .cvgui_window
  set ::cvgui::cv_values {}
  foreach c [run_cv list] {
    lappend ::cvgui::cv_values [format_value [run_cv colvar $c value] ]
  }
}

proc ::cvgui::update_frame { name molid op } {
  # name == vmd_frame
  # molid == molecule id of the newly changed frame
  # op == w
  global vmd_frame
  variable ::cvgui::plothandle

  if { $molid != [molinfo top] } {
    return
  }
  set f [molinfo $molid get frame]

  set ::cvgui::current_frame $f
  run_cv frame $f
  ::cvgui::refresh_values

  if [info exists plothandle] {
    # detect if plot was closed
    if [catch {$plothandle getpath}] {
      unset plothandle
    } else {
      # delete previous lines by tinkering with Multiplot's internals
      set ns [namespace qualifiers $::cvgui::plothandle]
      set ${ns}::vline {}
      $plothandle configure -vline [list $f -dash "-"] -plot
    }
  }
}


proc ::cvgui::format_value val {
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


proc ::cvgui::selected {} {
  set w .cvgui_window
  set l {}
  refresh_table ;# to make sure selected colvars haven't been deleted
  foreach i [$w.cvtable curselection] {
    lappend l [lindex $::cvgui::cvs $i]
  }
  return $l
}


# Enable or disable real-time tracking of VMD frame
proc ::cvgui::change_track_frame {} {
  global vmd_frame

  set molid [molinfo top]
  if {$::cvgui::track_frame} {
    trace add variable vmd_frame($molid) write ::cvgui::update_frame
    update_frame vmd_frame [molinfo top] w
  } else {
    trace remove variable vmd_frame($molid) write ::cvgui::update_frame
  }
}

# Load config from file
proc ::cvgui::load {} {
  set path [tk_getOpenFile -filetypes {{"Colvars cfg" .in} {"Colvars cfg" .colvars}}]
  if [string compare $path ""] {
    run_cv configfile $path
    refresh_table
  }
}


# Save config of colvars to file (can we do it w/ biases ? need bias type keyword)
proc ::cvgui::save {} {

  set path [tk_getSaveFile -filetypes {{"Colvars cfg" .in} {"Colvars cfg" .colvars}}]

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


proc ::cvgui::del {} {
  foreach c [selected] {
    run_cv colvar $c delete
  }
  refresh_table
}


proc ::cvgui::add {} {
  edit true
}


# Enable or disable real-time tracking of VMD frame
proc ::cvgui::edit { {add false} } {
  set cfg ""

  if $add {
    # do not remove existing vars
    set cvs {}
    set cfg "colvar {\n  name <x>\n  <component> {\n    <atomGroup> { }\n  }\n}\n"
  } else {
    set cvs [selected]
    foreach c $cvs {
      append cfg "colvar {"
      append cfg [run_cv colvar $c getconfig]
      append cfg "}\n\n"
    }
    set ::cvgui::backup_cfg $cfg
  }
  set w .cvgui_window
  set editor [toplevel $w.editor]
  tk::text $w.editor.text
  $w.editor.text insert 1.0 $cfg
  set ::cvgui::being_edited $cvs

  grid $w.editor.text -row 0 -columnspan 2 -sticky nsew

  ttk::button $w.editor.apply -text "Apply" -command ::cvgui::edit_apply
  ttk::button $w.editor.cancel -text "Cancel" -command ::cvgui::edit_cancel
  grid $w.editor.apply -row 1 -column 0 -sticky e
  grid $w.editor.cancel -row 1 -column 1 -sticky w
  grid columnconfigure $w.editor 0 -weight 1
  grid columnconfigure $w.editor 1 -weight 1
  grid rowconfigure $w.editor 0 -weight 1
}


proc ::cvgui::edit_apply {} {
  set w .cvgui_window
  foreach c $::cvgui::being_edited {
    puts "Deleting $c"
    run_cv colvar $c delete
    puts [cv list]
  }
  set cfg [$w.editor.text get 1.0 end-1c]
  if { $cfg != "" } {
    set res [run_cv config $cfg]
    if { [string compare $res ""] } {
      # error: restore backed up cfg
      run_cv config $::cvgui::backup_cfg
    }
  }
  set ::cvgui::being_edited {}
  destroy $w.editor
  refresh_table
}


proc ::cvgui::reset {} {
  run_cv reset
  refresh_table
}


proc ::cvgui::edit_cancel {} {
  set w .cvgui_window
  set ::cvgui::being_edited {}
  destroy $w.editor
}


proc ::cvgui::plot {} {
  variable ::cvgui::plothandle
  set nf [molinfo top get numframes]
  set x {}
  for {set f 0} {$f < $nf} {incr f} { lappend x $f }

  set cvs [selected]
  if { [llength $cvs] == 0 } {
    # If no selection, plot all variables
    set cvs [run_cv list]
    if { [llength $cvs] == 0 } {
      return
    }
  }

  foreach c $cvs {
    set val [run_cv colvar $c update]
    set size [llength $val]
    if { $size == 1 } {
      set names($c) $c
      set y($c) {}
    } else {
      for {set i 1} {$i <= $size} {incr i} {
        set n $c; append n "_$i"
        lappend names($c) $n
        set y($n) {}
      }
    }
  }

  for {set f 0} {$f< $nf} {incr f} {
    run_cv frame $f
    foreach c $cvs {
      set val [run_cv colvar $c update]
      foreach ni $names($c) vi $val {
        lappend y($ni) $vi
      }
    }
  }
  set plothandle [multiplot -title "Colvars trajectory" -xlabel "Frame" -ylabel "Colvar value"]
  foreach c $cvs {
    foreach n $names($c) {
      $plothandle add $x $y($n) -legend $n
    }
  }
  $plothandle configure -vline [list $::cvgui::current_frame -dash "-"] -plot
  $plothandle replot
}


# The main window
proc ::cvgui::createWindow {} {
  set w [toplevel .cvgui_window]

  # setup Colvars if not already there
  if [catch { cv version}] {
    run_cv molid top
  }
  set gridrow 0
  grid [ttk::button $w.load -text "Load file" -command ::cvgui::load -padding "2 0 2 0"] -row $gridrow -column 0 -pady 5 -padx 2 -sticky nsew
  grid [ttk::button $w.save -text "Save colvars" -command ::cvgui::save -padding "2 0 2 0"] -row $gridrow -column 1 -pady 5 -padx 2 -sticky nsew
  grid [ttk::button $w.reset -text "Reset" -command ::cvgui::reset -padding "2 0 2 0"] -row $gridrow -column 2 -pady 5 -padx 2 -sticky nsew

  tablelist::tablelist $w.cvtable -columns {
    0 Colvars
  } -stretch all -selectmode extended -selecttype row -listvariable ::cvgui::cvs
  tablelist::tablelist $w.valuetable -columns {
    0 Values
  } -stretch all -listvariable ::cvgui::cv_values

  refresh_table

  incr gridrow
  grid $w.cvtable -row $gridrow -column 0 -sticky news
  grid $w.valuetable -row $gridrow -column 1 -sticky news -columnspan 2
  grid rowconfigure $w $gridrow -weight 1

  incr gridrow
  grid [ttk::button $w.edit -text "Edit" -command ::cvgui::edit -padding "2 0 2 0"] -row $gridrow -column 0 -pady 5 -padx 2 -sticky nsew
  grid [ttk::button $w.add -text "Add" -command ::cvgui::add -padding "2 0 2 0"] -row $gridrow -column 1 -pady 5 -padx 2 -sticky nsew
  grid [ttk::button $w.del -text "Delete" -command ::cvgui::del -padding "2 0 2 0"] -row $gridrow -column 2 -pady 5 -padx 2 -sticky nsew

  incr gridrow
  grid [ttk::button $w.refresh -text "Refresh" -command ::cvgui::refresh_table -padding "2 0 2 0"] -row $gridrow -column 0 -pady 5 -padx 2 -sticky nsew
  grid [ttk::button $w.plot -text "Plot" -command ::cvgui::plot -padding "2 0 2 0"] -row $gridrow -column 1 -pady 5 -padx 2 -sticky nsew
  #grid [ttk::button $w.find -text "Find frame" -command ::cvgui::find -padding "2 0 2 0"] -row $gridrow -column 2 -pady 5 -padx 2 -sticky nsew

  incr gridrow
  grid [label $w.frameTxt -text "Frame:"] -row $gridrow -column 0 -pady 5 -padx 2 -sticky nsew
  grid [label $w.frame -textvariable ::cvgui::current_frame] -row $gridrow -column 1 -pady 5 -padx 2 -sticky nsew
  grid [ttk::checkbutton $w.trackFrame -text "Track" -command ::cvgui::change_track_frame -variable ::cvgui::track_frame]  -row $gridrow -column 2  -pady 5 -padx 2 -sticky nsew
  change_track_frame ;# activate tracking if necessary

  grid columnconfigure $w 0 -weight 1
  grid columnconfigure $w 1 -weight 1
  grid columnconfigure $w 2 -weight 1
}

::cvgui::createWindow
