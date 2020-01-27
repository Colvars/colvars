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
# which is most consistent for trajectory animation (determined by the frame number of mol)

# TODO Multiplot:
# - properly calculate position of cursor in plot when not all the plot is visible (resized window)
# - handle several windows at once - at least one pairwise, one timeline
# - integrate interactive hacks into interface
# - display pairwise traj on top of known 2D data (eg. FE surface)

# TODO maybe:
# - index group builder

package provide cv_dashboard 1.3

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

  variable atom_rep     ;# hash array of: list of macro names, list of atom representations, indexed by colvar name
  variable grad_objects ;# hash array ids of graphical objects displaying gradients, indexed by colvar name
  variable grad_scale_choice ;# whether to scale gradients by a fixed factor of to obtain given max norm
  variable grad_norm 5.0  ;# Default value for gradient max norm
  variable grad_scale 1.0 ;# Default value for gradient scale factor

  variable mol -1       ;# ID of molecule currently associated with Colvars

  variable indent "    " ;# indentation for config strings

  variable units
  variable units_to_text
  variable text_to_units
  array set text_to_units {"real (Angstrom, kcal/mol)" "real" "Gromacs (nm, kJ/mol)" "gromacs" \
    "metal (Angstrom, eV)" "metal" "electron (Bohr, Hartree)" "electron"}
  # Build reverse map
  foreach { text units } [array get text_to_units] {
    set units_to_text($units) $text
  }

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
  if {[molinfo num] == 0 } {
    tk_messageBox -icon error -title "Colvars Dashboard Error"\
      -message "No molecule loaded. Please load a molecule and use the Reset button.\n"
  } else {

    if [catch { cv version }] {
      # if not already there, setup Colvars on top molecule
      ::cv_dashboard::run_cv molid [molinfo top]
    }
    # Try to obtain actual molid, if that fails default to top
    if [catch { set ::cv_dashboard::mol [cv molid]} ] {
      set ::cv_dashboard::mol [molinfo top]]
    }
  }

  if {[winfo exists .cv_dashboard_window]} {
    wm deiconify .cv_dashboard_window
    return .cv_dashboard_window
  }

  return [eval ::cv_dashboard::createWindow]
}


#################################################################
# General utilities, Colvars Module related
#################################################################


# Call the "cv" interface to Colvars, catching errors and displaying them to the user
proc ::cv_dashboard::run_cv args  {
  if { [lindex $args 0] != "molid" } {
    # Try to initialize the module if not there yet
    if [catch {cv version}] {
      if [catch {cv molid $::cv_dashboard::mol}] {
        # If that didn't work, don't try to proceed
        return
      }
    }
  }
  if [catch { cv {*}$args } res] {
    set short_cmd [string range $args 0 200]
    set short_message [string range $res 0 200]
    tk_messageBox -icon error -title "Colvars error" -parent .cv_dashboard_window\
      -message "Error running command:\n$short_cmd" -detail "$short_message\n\nSee console for further details."
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

  set cfg [substitute_atomselects $cfg]

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

  # Update atom visualizations for modified colvars
  foreach cv [dict keys $cv_configs] {
    if { [info exists ::cv_dashboard::atom_rep($cv)] } {
      show_atoms $cv
    }
  }

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
  refresh_units
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
      if { [regexp -nocase {^\s*colvar\s+\{\s*(.*)} $line match firstline] } {
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
    # non-word characters are spaces and {}# (do not use Tcl's restrictive \w)
    if { $brace_depth == 1 } {
      regexp -nocase {^\s*name\s+([^\s{}#]+)} $line match name
    }

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


# Detect magic comments and update following atom keywords
proc ::cv_dashboard::substitute_atomselects { cfg_in } {
  set lines [split $cfg_in "\n"]
  set cfg ""
  set seltext ""
  foreach line $lines {
    # If we:
    # 1) remember a seltext from previous line, and
    # 2) are processing a line starting with "atom" (atom selection keyword)
    if { $seltext != "" && [regexp -nocase {^(\s*)atom} $line match spaces] } {
      set sel [atomselect $::cv_dashboard::mol $seltext]
      set serials [$sel get serial]
      $sel delete
      if {[llength $serials] == 0 } {
        tk_messageBox -icon error -title "Colvars warning" -parent .cv_dashboard_window\
          -message "Selection text \"${seltext}\" for automatic atom selection matches zero atoms. \
Keeping atom numbers from existing configuration."
        # Keep existing atom definition line
        append cfg $line "\n"
        # Forget seltext for next lines
        set seltext ""
      } else {
        # Replace keyword, keeping indenting spaces
        append cfg "${spaces}atomNumbers $serials\n"
        # Forget seltext for next lines
        set seltext ""
      }
    } else {
      append cfg $line "\n"
    }
    regexp {^\s*# auto-updating selection: "(.*)"} $line match seltext
  }
  return $cfg
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


# Callback to update CV Dashboard when VMD's molecule changes to new frame
proc ::cv_dashboard::update_frame { name molid op } {
  # name == vmd_frame
  # molid == molecule id of the newly changed frame
  # op == w

  if { $molid != $::cv_dashboard::mol } {
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


# React to molecules being created or deleted
proc ::cv_dashboard::update_mol_list { name molid op } {
  # Update mol indices in combobox
  .cv_dashboard_window.mol configure -values [molinfo list]

  # Did we just lose the molecule Colvars was connected to?
  if { ($molid == $::cv_dashboard::mol) && ($::vmd_initialize_structure($molid) == 0) } {
    tk_messageBox -icon error -title "Colvars Dashboard Error"\
      -message "The molecule associated to the Colvars Module was deleted.\nSave the configuration if necessary, load a molecule and use the Reset button.\n"
    # remove tracking of deleted molecule
    trace remove variable ::vmd_frame($molid) write ::cv_dashboard::update_frame
    set ::cv_dashboard::mol -1
  }

  # Did we just add a molecule while we had none available? Default to that
  if { ($::cv_dashboard::mol == -1) && ($::vmd_initialize_structure($molid) == 1)} {
    set ::cv_dashboard::mol $molid
    .cv_dashboard_window.mol set $molid
  }
}


proc ::cv_dashboard::change_mol {} {
  set newmolid [.cv_dashboard_window.mol get]

  if { $newmolid != $::cv_dashboard::mol } {
    trace remove variable ::vmd_frame($::cv_dashboard::mol) write ::cv_dashboard::update_frame
    # Remove all graphical objects which would be orphaned
    ::cv_dashboard::hide_all_atoms
    ::cv_dashboard::hide_all_gradients

    set ::cv_dashboard::mol $newmolid
    # Remember config
    if {$::cv_dashboard::units == ""} {
      set cfg ""
    } else {
      set cfg "units $::cv_dashboard::units\n\n"
    }
    foreach cv [run_cv list] {
        append cfg "colvar {[get_config $cv]}\n\n"
    }
    reset
    apply_config $cfg
    change_track_frame ;# activate tracking of new molecule if requested
  }
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
