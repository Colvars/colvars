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

package provide cv_dashboard 1.1

namespace eval ::cv_dashboard {
  # General UI state
  variable current_frame 0  ;# linked to frame display
  variable cvs {}           ;# linked to colvar table
  variable track_frame 1    ;# start tracking by default

  # State variables for config editor
  variable being_edited
  variable backup_cfg
  variable filetype "atomsFile"
  variable colvar_configs  ;# dictionary mapping names to cfg strings
  set colvar_configs [dict create]

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
    set script_dir [file dirname [info script]]
    set template_dir ${script_dir}/templates
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
# General utilities, Colvars Module related
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


# Apply a new config string:
# - dump a file with debug information to help if an wrong config causes a crash
# - save list of existing colvars
# - submit to Colvars module for parsing
# - extract config strings of individual colvars
# - keep full config string with comments for any newly added colvar
proc ::cv_dashboard::apply_config { cfg } {
  if { $cfg == "" } {
    return ""
  }
  # Dump config for debugging possible crashes
  set dump [open "_dashboard_saved_config.colvars" w]
  puts $dump "# Current configuration of Colvars Module\n"
  foreach c [run_cv list] {
      puts $dump "colvar {[get_config $c]}\n"
  }
  puts $dump "\n# New config string to be applied\n"
  puts $dump $cfg
  close $dump

  set cvs_before [run_cv list]
  # Actually submit new config to the Colvars Module
  set res [run_cv config $cfg]
  set cvs_after [run_cv list]

  # Extract config for individual colvars
  set cv_configs [extract_colvar_configs $cfg]

  # Completely update the map of colvar configs
  set new_map [dict create]
  dict for { name cfg } $cv_configs {
    # Only record config for cvs that actually appeared just now
    if { ([lsearch $cvs_after $name] > -1) && ([lsearch $cvs_before $name] == -1) } {
      dict set new_map $name $cfg
    }
  }
  # Look for missing cvs in the old map
  foreach cv $cvs_after {
    if { ! [dict exists $new_map $cv]} {
      catch {
        dict set new_map $cv [dict get $::cv_dashboard::colvar_configs $cv]
      }
    }
  }
  # Overwrite old map
  set ::cv_dashboard::colvar_configs $new_map
  refresh_table
  return $res
}


# Parse config string to extract colvar blocks
# Return dictionary of colvar names -> config strings
# Needs to fail gracefully upon unmatched braces
proc ::cv_dashboard::extract_colvar_configs { cfg_in } {
  set lines [split $cfg_in "\n"]
  set in_cv 0
  set brace_depth 0
  set map [dict create]
  set name ""
  foreach line $lines {
    if { $in_cv == 0 } {
      # In main body, just look for colvar definition
      if { [regexp {^\s*colvar\s+\{\s*(.*)} $line match firstline] } {
        set in_cv 1
        set cv_line 1
        set brace_depth 1
        set cv_cfg "\n"
        # The first line may follow the opening brace immediately
        if { [string length $firstline] } {
          set line "    ${firstline}"
        } else {
          # Nothing more to parse from this line
          continue
        }
      } else {
        # Don't parse non-colvar data
        continue
      }
    }
    # Now we're parsing a line of colvar config, try to get name
    regexp {^\s*name\s+([^\s]+)} $line match name

    # Finally, the tedious fishing for braces
    regexp {^[^#]*} $line nocomments
    set chars [split $nocomments ""]
    set cur_line ""
    foreach c $chars {
      switch $c {
        "{" { incr brace_depth }
        "}" {
          incr brace_depth -1
          if { $brace_depth < 0 } {
            # probably mismatched braces
            # give up on parsing the rest but try to return any variable already parsed
            puts "Warning: mismatched braces in configuration line:\n${line}"
            return $map
          }
          if { $brace_depth == 0 } {
            # End of colvar block
            if { [string length $cur_line] > 0 } {
              append cv_cfg $cur_line "\n"
            }
            dict set map $name $cv_cfg
            set in_cv 0
            set name ""
          }
        }
      }
      # keep track of line up to current char to catch the last line
      append cur_line $c
    }

    if { $in_cv } {
      if { $cv_line >= 1 } {
        append cv_cfg $line "\n"
      }
      incr cv_line
    }
  }
  return $map
}


# Looks for config in saved map
# if not found, queries the colvar itself (get stripped cfg)
proc ::cv_dashboard::get_config { cv } {
  if { [dict exists $::cv_dashboard::colvar_configs $cv] } {
    return [dict get $::cv_dashboard::colvar_configs $cv]
  } else {
    return [run_cv colvar $cv getconfig]
  }
}


#################################################################
# GUI-related utilities
#################################################################


# Callback to update CV Dashboard when VMD's top molecule changes to new frame
proc ::cv_dashboard::update_frame { name molid op } {
  # name == vmd_frame
  # molid == molecule id of the newly changed frame
  # op == w

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
  set h [toplevel $parent.helpWindow]
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


# Add keyboard bindings for trajectory animation to widget given by path
proc ::cv_dashboard::traj_animation_bindings { path } {
  bind $path <Left>           { ::cv_dashboard::chg_frame -1 }
  bind $path <Right>          { ::cv_dashboard::chg_frame 1 }
  bind $path <Shift-Left>     { ::cv_dashboard::chg_frame -10 }
  bind $path <Shift-Right>    { ::cv_dashboard::chg_frame 10 }
  bind $path <Control-Left>   { ::cv_dashboard::chg_frame -50 }
  bind $path <Control-Right>  { ::cv_dashboard::chg_frame 50 }
  bind $path <Home>           { ::cv_dashboard::chg_frame start }
  bind $path <End>            { ::cv_dashboard::chg_frame end }
}
