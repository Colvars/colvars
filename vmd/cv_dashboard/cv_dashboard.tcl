# Colvars Dashboard -- based on the Colvars Module for VMD
# Jérôme Hénin <henin@ibpc.fr> 2018

# Usage (after installing):
# package require cv_dashboard
# cv_dashboard

# Design principles:
# - take advantage of colvars/VMD binding for maximum user interaction
# - hide the colvars config text from user, instead expose colvar, names and values
# - do not try to parse the colvars config (let the Colvars Module do it)
#   to avoid coming up with an incompatible parser

# This plugin only acts on the "top" molecule
# which is most consistent for trajectory animation (determined by the frame number of top mol)

# TODO Multiplot:
# - properly calculate position of cursor in plot when not all the plot is visible (resized window)
# - handle several windows at once - at least one pairwise, one timeline
# - integrate interactive hacks into interface
# - display pairwise traj on top of known 2D data (eg. FE surface)

# TODO maybe:
# - index group builder

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

set script_dir [file dirname [info script]]
source [file join $script_dir cv_dashboard_main.tcl]
source [file join $script_dir cv_dashboard_editor.tcl]
source [file join $script_dir cv_dashboard_plots.tcl]


proc cv_dashboard {} {
  return [eval ::cv_dashboard::createWindow]
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


# Displays a non-blocking help window with the provided info
proc ::cv_dashboard::help_window { parent wtitle title text } {
  set h [toplevel $parent.help]
  wm title $h $wtitle
  tk::text $h.text -yscrollcommand [list $h.vsb set]
  ttk::scrollbar $h.vsb -orient vertical -command [list $h.text yview]

  $h.text insert insert ${title}\n\n title
  $h.text tag configure title -font "Helvetica -14 bold" -justify center
  $h.text insert insert $text
  $h.text configure -state disabled
  ttk::button $h.close -text "Close" -command " destroy $h " -padding "2 0 2 0"

  grid $h.text -row 0 -column 0 -sticky nsew
  grid $h.vsb -row 0 -column 1 -sticky nsew
  grid $h.close -row 1 
  grid columnconfigure $h 0 -weight 1
  grid rowconfigure $h 0 -weight 1
}
