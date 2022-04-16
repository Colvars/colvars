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

package provide cv_dashboard 1.5

namespace eval ::cv_dashboard {
  # General UI state
  variable current_frame 0  ;# linked to frame display
  variable cvs {}           ;# linked to colvar table
  variable track_frame 1    ;# start tracking by default

  # State variables for config editor
  variable being_edited
  variable backup_cfg
  variable filetype        "atomsFile"
  variable colvar_configs  [dict create] ;# dictionary mapping names to cfg strings
  variable bias_configs    [dict create]
  variable global_config   [dict create]
  variable global_comments ""
  variable cv_clipboard    [dict create]
  variable bias_clipboard    [dict create]

  # Handle to keep track of interactive plot
  variable plothandle
  variable plottype      ;# timeline, 2cv, histogram

  variable atom_rep      ;# hash array of: list of macro names, list of atom representations, indexed by colvar name
  variable grad_objects  ;# hash array ids of graphical objects displaying gradients, indexed by colvar name
  variable force_objects ;# hash array ids of graphical objects displaying forces, indexed by bias name
  variable grad_scale_choice ;# whether to scale gradients by a fixed factor of to obtain given max norm
  variable grad_norm 5.0  ;# Default value for gradient max norm
  variable grad_scale 1.0 ;# Default value for gradient scale factor

  variable mol -1       ;# ID of molecule currently associated with Colvars

  variable indent "    " ;# indentation for config strings
  variable font   {normal 9} ;# default font for labels

  variable units
  variable units_to_text
  variable text_to_units
  array set text_to_units {"real (Angstrom, kcal/mol)" "real" "Gromacs (nm, kJ/mol)" "gromacs" \
    "metal (Angstrom, eV)" "metal" "electron (Bohr, Hartree)" "electron"}
  # Build reverse map
  foreach { text units } [array get text_to_units] {
    set units_to_text($units) $text
  }
  dict set global_config units "real"

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
source [file join $script_dir cv_dashboard_display.tcl]
source [file join $script_dir cv_dashboard_settings.tcl]


proc cv_dashboard {} {

  if { [vmdinfo versionmsg] == "VMD for LINUXAMD64, version 1.9.3 (November 30, 2016)" } {
    tk_messageBox -icon error -title "Colvars Dashboard Error"\
      -message "VMD 1.9.3 (released 2016) ships with an unmaintained version of the Colvars Module.\n
Please upgrade to VMD 1.9.4 alpha or later."
    return
  }

  if {[molinfo num] == 0 } {
    tk_messageBox -icon error -title "Colvars Dashboard Error"\
      -message "No molecule loaded. Please load a molecule and use the Reset button.\n"
  } else {

    set cv_mol -1
    catch { set cv_mol [cv molid] }
    if { $cv_mol == -1 } {
      # if not already there, setup Colvars on top molecule
      ::cv_dashboard::run_cv molid [molinfo top]
    }
    # Try to obtain actual molid, if that fails (older cv API) default to top
    if [catch { set ::cv_dashboard::mol [cv molid]} ] {
      set ::cv_dashboard::mol [molinfo top]
    }
  }

  if {[winfo exists .cv_dashboard_window]} {
    wm deiconify .cv_dashboard_window
    return .cv_dashboard_window
  }

  if {[llength [info procs calc_colvar_forces]] == 0} {
    # Create dummy proc to avoid error messages is scriptedColvarForces is enabled
    proc calc_colvar_forces { ts } {}
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

  # Don't try to update if no frames are loaded
  if { $args == "update" } {
    set nf [molinfo $::cv_dashboard::mol get numframes]
    if { $nf <= 0 } {
      puts "No frames loaded, cannot update Colvars"
      return
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

  # Save viewpoints for all molecules, as they could be reset when applying config
  set vp [get_viewpoints]

  set cvs_before [run_cv list]
  set biases_before [run_cv list biases]
  # Actually submit new config to the Colvars Module
  set res [run_cv config $cfg]
  set cvs_after [run_cv list]
  set biases_after [run_cv list biases]

  set_viewpoints $vp

  # Extract config for individual colvars and biases
  lassign [extract_configs $cfg] cv_configs bias_configs global_config comments

  if { $comments != "" } {
    append ::cv_dashboard::global_comments "${comments}\n"
  }

  set ::cv_dashboard::global_config [dict merge $::cv_dashboard::global_config $global_config]

  # Update atom visualizations for modified colvars
  foreach cv [dict keys $cv_configs] {
    if { [info exists ::cv_dashboard::atom_rep($cv)] } {
      show_atoms $cv
    }
  }

  # Update the map of colvar configs
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
  # Overwrite old colvar map
  set ::cv_dashboard::colvar_configs $new_map

  # Update the map of bias configs
  set new_bias_map [dict create]
  dict for { name cfg } $bias_configs {
    # Only record config for biases that actually appeared just now
    if { ([lsearch $biases_after $name] > -1) && ([lsearch $biases_before $name] == -1) } {
      dict set new_bias_map $name $cfg
    }
  }
  # Look for missing biases in the old map
  foreach bias $biases_after {
    if { ! [dict exists $new_bias_map $bias]} {
      catch {
        dict set new_bias_map $bias [dict get $::cv_dashboard::bias_configs $bias]
      }
    }
  }
  # Overwrite old bias map
  set ::cv_dashboard::bias_configs $new_bias_map

  refresh_table
  refresh_units
  # Refresh map in case list of Colvars atoms has changed
  catch { unset ::cv_dashboard::atom_id_map }
  return $res
}

# Parse config string to extract colvar blocks
# Return dictionary of colvar names -> config strings
# Needs to fail gracefully upon unmatched braces
proc ::cv_dashboard::extract_configs { cfg_in } {

  set indent $::cv_dashboard::indent

  set biases [cv list biases]
  # "cv bias name type" after 2021-12-07
  if { [string compare [run_cv version] "2021-12-20"] >= 0 } {
    set bias_types [list]
    foreach cvb $biases {
      lappend bias_types [string tolower [cv bias $cvb type]]
    }
    set bias_types [lsort -unique $bias_types]
  } else {
    # if not available, use hard-coded list
    set bias_types [list abf alb harmonic harmonicwalls histogram histogramrestraint linear metadynamics reweightamd]
  }
  foreach t $bias_types {
    set anonymous_bias_count($t) 0
  }

  set lines [split $cfg_in "\n"]
  set in_block 0
  set brace_depth 0
  set cv_map [dict create]        ;# cv name -> config (with comments)
  set bias_map [dict create]      ;# bias name -> config (with comments)
  set global_cfg_map [dict create] ;# keyword -> rest of the line (value + comments)
  set comment_lines ""            ;# lines with only comments
  set name ""
  set keyword ""

  foreach line $lines {
    if { $in_block == 0 } {
      # In main body, look for block definition
      if { [regexp -nocase {^\s*(\S+)\s+\{\s*(.*)} $line match keyword firstline] } {
        set keyword [string tolower $keyword]
        set in_block 1
        set block_line 1
        set brace_depth 1
        set block_cfg "\n"
        # The first line may follow the opening brace immediately
        if { [string length $firstline] } {
          set line "${indent}${firstline}"
        } else {
          # Nothing more to parse from this line
          continue
        }
      } else {
        # Not a block; it's either a comment line...
        if {[regexp -nocase {^\s*\#} $line] } {
          append comment_lines "${line}\n"
        } elseif { [regexp -nocase {^\s*(\S+)\s+(.*)} $line match keyword value] } {
          # or it goes to general Colvars module config
          dict set global_cfg_map $keyword $value
        }
        continue
      }
    }
    # Now we're parsing a line of block (colvar/bias) config, try to get name
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
            return $cv_map
          }
          if { $brace_depth == 0 } {
            # End of block
            if { [string length $cur_line] > 0 } {
              append block_cfg $cur_line "\n"
            }
            if {$keyword == "colvar"} {
              dict set cv_map $name $block_cfg
            } else {
              if { [lsearch $bias_types $keyword] > -1 } {
                # Bias names are unique, but not always available from the config
                if { $name == "" } {
                  incr anonymous_bias_count($keyword)
                }
                dict set bias_map [list $keyword $anonymous_bias_count($keyword) $name] $block_cfg
              } else {
                # What to make of unrecognized block keyword?
                puts "Warning: unrecognized block keyword in input: $keyword"
                puts "with configuration block:"
                puts "---------\n$block_cfg\n---------\n"
              }
            }
            set in_block 0
            set name ""
            set keyword ""
          }
        }
      }
      # keep track of line up to current char to catch the last line
      append cur_line $c
    }

    if { $in_block } {
      if { $block_line >= 1 } {
        append block_cfg $line "\n"
      }
      incr block_line
    }
  }

  set new_bias_map [dict create]
  # Done parsing, now find out missing bias names
  foreach key [dict keys $bias_map] {
    lassign $key keyword index name
    if { $name == "" } {
      # find generated bias name
      # the last n biases with names "keyword$i" are the ones we want
      set auto_names [regexp -all -inline -nocase "$keyword\\d+" $biases]
      set id [expr {[llength $auto_names] -1 - $anonymous_bias_count($keyword) + $index}]
      set name [lindex $auto_names $id]
    }
    # New dict has name as key and "keyword { cfg } " as value
    dict set new_bias_map $name [list $keyword [dict get $bias_map $key]]
  }
  return [list $cv_map $new_bias_map $global_cfg_map $comment_lines]
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
      if {[$sel num] == 0 } {
        tk_messageBox -icon error -title "Colvars warning" -parent .cv_dashboard_window\
          -message "Selection text \"${seltext}\" for automatic atom selection matches zero atoms. \
Keeping atom numbers from existing configuration."
        # Keep existing atom definition line
        append cfg $line "\n"
        # Forget seltext for next lines
        set seltext ""
      } else {
        # Replace keyword, keeping indenting spaces
        append cfg "${spaces}[sel2cvatoms $sel]\n"
        # Forget seltext for next lines
        set seltext ""
      }
      $sel delete
    } else {
      append cfg $line "\n"
    }
    regexp {^\s*# auto-updating selection: "(.*)"} $line match seltext
  }
  return $cfg
}


# Looks for config in saved map
# if not found, queries the colvar itself (get stripped cfg)
proc ::cv_dashboard::get_cv_config { cv } {
  if { [dict exists $::cv_dashboard::colvar_configs $cv] } {
    return [dict get $::cv_dashboard::colvar_configs $cv]
  } else {
    return [run_cv colvar $cv getconfig]
  }
}

# Looks for config in saved map
# if not found, queries the bias itself (get stripped cfg)
proc ::cv_dashboard::get_bias_keyword_config { bias } {
  # First try: look up in table of configs
  if { [dict exists $::cv_dashboard::bias_configs $bias] } {
    set key_cfg [dict get $::cv_dashboard::bias_configs $bias]
    # Cfg already includes bias keyword
    return $key_cfg
  # Second chance: query type and config via script interface
  } elseif { [string compare [run_cv version] "2021-12-20"] >= 0 } {
    set cfg [run_cv bias $bias getconfig]
    return [list [run_cv bias $bias type] $cfg]
  } else {
    puts "Warning: could not find configuration for bias $bias.\n"
  }
}

# Returns reconstructed config file for colvars, biases, the module, + comments
proc ::cv_dashboard::get_whole_config { } {

  set cfg $::cv_dashboard::global_comments
  if { $cfg != "" } {append cfg "\n"}

  dict for {key value} $::cv_dashboard::global_config {
    # Exclude index files, written separately below
    # Eack keyword is kept only once
    if { [string tolower $key] != "indexfile" } {
      append cfg "$key $value\n"
    }
  }

  set indexFiles [list]
  catch { set indexFiles [cv listindexfiles] }
  foreach ndx $indexFiles {
    append cfg "indexFile $ndx\n"
  }

  if { $cfg != "" } {append cfg "\n"}

  foreach cv [run_cv list] {
    append cfg "colvar {[get_cv_config $cv]}\n\n"
  }

  foreach bias [run_cv list biases] {
    lassign [get_bias_keyword_config $bias] keyword config

    # Skip if bias config was not found
    if { $keyword == {} } { continue }

    if { $keyword == "harmonicwalls" } {
      regexp -line -nocase {^\s*colvars\s+(.*)} $config match cvs
      if { $bias == "${cvs}w" } {
        continue
      }
    }
    append cfg "$keyword {$config}\n\n"
  }
  return $cfg
}

# Checks whether cv is associated to a volmap
proc ::cv_dashboard::is_volmap { cv } {
  catch {set id [cv colvar $cv getvolmapids]}
  if { [llength $id] > 1 } { return 1 }
  return [expr [lindex $id 0] != -1]
}


# Checks whether cv is a unit quaternion (current value)
proc ::cv_dashboard::is_unit_quaternion { cv } {
  set EPSILON 1e-10
  set v [run_cv colvar $cv value]
  if {[llength $v] != 4} { return 0 }
  set l [expr abs([veclength $v] - 1.0)]
  return [expr $l < $EPSILON]
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

  if { $f < 0 } { return }
  # set Colvars Module to requested frame
  run_cv frame $f
  # refresh dashboard table
  refresh_values
  # refresh the frame marker in the plot
  display_marker $f
  # refresh displayed CV gradients and bias forces
  update_shown_gradients
  update_shown_forces
  # refresh any displayed rotation operators
  update_rotation_display
}


# React to molecules being created or deleted
proc ::cv_dashboard::update_mol_list { name molid op } {
  set main .cv_dashboard_window.tabs.main

  # Update mol indices in combobox
  $main.mol configure -values [molinfo list]

  # Did we just lose the molecule Colvars was connected to?
  if { ($molid == $::cv_dashboard::mol) && ($::vmd_initialize_structure($molid) == 0) } {
    tk_messageBox -icon error -title "Colvars Dashboard Error"\
      -message "The molecule associated to the Colvars Module was deleted.\nSave the configuration if necessary before loading a molecule.\n"
    # remove tracking of deleted molecule
    trace remove variable ::vmd_frame($molid) write ::cv_dashboard::update_frame
    set ::cv_dashboard::mol -1
  }

  # Did we just add a molecule while we had none available? Default to that
  if { ($::cv_dashboard::mol == -1) && ($::vmd_initialize_structure($molid) == 1)} {
    set ::cv_dashboard::mol $molid
    $main.mol set $molid
    catch {cv delete}
    run_cv molid $molid
    change_track_frame ;# activate tracking of new molecule if requested
    refresh_table
  }
}


proc ::cv_dashboard::change_mol {} {
  set main .cv_dashboard_window.tabs.main
  set newmolid [$main.mol get]

  if { $newmolid != $::cv_dashboard::mol } {
    trace remove variable ::vmd_frame($::cv_dashboard::mol) write ::cv_dashboard::update_frame
    # Remove all graphical objects which would be orphaned
    ::cv_dashboard::hide_all_atoms
    # Backup shown grads/forces
    set grad_objects [array get ::cv_dashboard::grad_objects]
    set force_objects [array get ::cv_dashboard::force_objects]
    ::cv_dashboard::hide_all_gradients
    ::cv_dashboard::hide_all_forces

    set ::cv_dashboard::mol $newmolid

    # Remember config
    set cfg [get_whole_config]

    reset
    apply_config $cfg
    # Attempt to restore grads/forces
    array set ::cv_dashboard::grad_objects $grad_objects
    array set ::cv_dashboard::force_objects $force_objects
    change_track_frame ;# activate tracking of new molecule if requested
    refresh_table
  }
}


proc ::cv_dashboard::switch_to_top_mol {} {
  set main .cv_dashboard_window.tabs.main
  $main.mol set [molinfo top]
  change_mol
}


# Displays a non-blocking help window with the provided info
proc ::cv_dashboard::help_window { parent wtitle title text } {
  catch { destroy $parent.helpWindow }
  set h [toplevel $parent.helpWindow]
  wm title $h $wtitle
  tk::text $h.text -yscrollcommand [list $h.vsb set] -bg white
  ttk::scrollbar $h.vsb -orient vertical -command [list $h.text yview]
  $h.text configure -cursor arrow

  $h.text insert insert ${title}\n\n title
  $h.text tag configure title -font "Helvetica 12 bold" -justify center
  if { [encoding system] != "utf-8" } {
    # Necessary to get special characters on Windows
    set text [encoding convertfrom "utf-8" $text]
  }
  $h.text insert insert $text

  add_tag $h.text {https?://\S+} URL
  $h.text tag configure URL -font "Mono 9" -foreground blue -underline true
  $h.text tag bind URL <Button-1> "::cv_dashboard::clickLink $h.text %x %y"

  add_tag $h.text {^#.*$} subtitle
  $h.text tag configure subtitle -font "Helvetica 12 bold"

  $h.text configure -state disabled
  ttk::button $h.close -text "Close" -command "destroy $h" -padding "2 0 2 0"

  grid $h.text -row 0 -column 0 -sticky nsew
  grid $h.vsb -row 0 -column 1 -sticky nsew
  grid $h.close -row 1
  grid columnconfigure $h 0 -weight 1
  grid rowconfigure $h 0 -weight 1
}

proc ::cv_dashboard::add_tag { t re tag } {
  set cnt 0
  set cur 0.0
  while 1 {
    set cur [$t search -count length -regexp -- $re $cur end]
    if {$cur == ""} { break }
    $t tag add $tag $cur "$cur + $length char"
    set cur [$t index "$cur + $length char"]
    incr cnt
  }
}

proc ::cv_dashboard::clickLink { text xpos ypos } {
  set i [$text index @$xpos,$ypos]
  set range [$text tag prevrange URL $i]
  set url [eval $text get $range]
  invokeBrowser $url
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


# Round floating-point number to $n significant figures
# without going to string representation, unlike format

proc ::cv_dashboard::round {x n} {
  if { $x == 0. } { return 0. }
  # e = 10^p, where p is the "number of decimal places" to keep
  # which can be negative
  set e [expr {pow(10, $n - floor(log10(abs($x))) - 1)}]
  return [expr {round($x * $e) / double($e)}]
}


proc ::cv_dashboard::get_viewpoints {} {
  set vp [dict create]
  # get the current matrices
  foreach mol [molinfo list] {
    dict set vp $mol [list \
      [molinfo $mol get rotate_matrix] \
      [molinfo $mol get center_matrix] \
      [molinfo $mol get scale_matrix]  \
      [molinfo $mol get global_matrix]]
  }
  return $vp
}


proc ::cv_dashboard::set_viewpoints { vp } {
  foreach mol [molinfo list] {
    if [dict exists $vp $mol] {
      lassign [dict get $vp $mol] a b c d
      molinfo $mol set rotate_matrix   $a
      molinfo $mol set center_matrix   $b
      molinfo $mol set scale_matrix    $c
      molinfo $mol set global_matrix   $d
    }
  }
}


# Find closest value in sorted list
proc ::cv_dashboard::bisect { lst val } {
  set len [llength $lst]
  set start 0
  set end [expr $len - 1]

  if { $len < 2 || $val <= [lindex $lst 0] } { return 0 }
  if { $val >= [lindex $lst $end] } { return $end }

  set mid [expr $len / 2]
  while { $start != $end } {
    if { [expr {$val <= [lindex $lst $mid]}] } {
      set end $mid
    } else {
      set start [expr {$mid + 1}]
    }
    set mid [expr {($start + $end ) / 2}]
  }
  if { $end == 0 } { return 0 }
  set start [expr $start - 1]

  set left [expr abs([lindex $lst $start] - $val)]
  set right [expr abs([lindex $lst $end] - $val)]
  if { $left <= $right } {
    return $start
  } else {
    return $end
  }
}


# Extract templates from files into dictionaries
proc ::cv_dashboard::parse_templates {} {
  foreach d { colvar component other bias } {
    set ::cv_dashboard::templates_$d [dict create]

    set path [file join $::cv_dashboard::template_dir "$d.colvars"]

    # Single-file template DBs
    if [catch {set db_file [open $path]}] {
      puts "No file $d"
      continue
    }

    set name ""
    while { [gets $db_file line] >= 0 } {
      if { [regexp "^#_(.+)" $line match newname] } {
        if { $name != "" } {
          dict set ::cv_dashboard::templates_$d $name $cfg
        }
        set name $newname
        set cfg ""
      } else {
        append cfg "$line\n"
      }
    }
    # Last config
    if { $name != "" } {
      dict set ::cv_dashboard::templates_$d $name $cfg
    }
    close $db_file
  }
}


# Create a variant of a provided name that is not in the list of reserved identifiers
proc ::cv_dashboard::make_unique_name { name reserved } {
  if { [lsearch $reserved $name] == -1 } { return $name }

  if { ![regexp {(.*)~(\d+)$} $name match base num] } {
    set base $name
    set num 0
  }

  set newname "${base}~${num}"
  while { [lsearch $reserved $newname] > -1 } {
    incr num
    set newname "${base}~${num}"
  }
  return $newname
}