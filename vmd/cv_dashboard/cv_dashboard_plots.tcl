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

  set nf [molinfo $::cv_dashboard::mol get numframes]
  # Get list of values for all frames
  for {set f 0} {$f< $nf} {incr f} {
    run_cv frame $f
    run_cv update
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

  traj_animation_bindings [set ${plot_ns}::w]
  if { $type == "timeline" } {
    bind [set ${plot_ns}::w] <Button-1>       { ::cv_dashboard::plot_clicked %x %y }
    bind [set ${plot_ns}::w] <Up>             { ::cv_dashboard::zoom 0.25 }
    bind [set ${plot_ns}::w] <Down>           { ::cv_dashboard::zoom 4 }
    bind [set ${plot_ns}::w] <Shift-Up>       { ::cv_dashboard::zoom 0.0625 }
    bind [set ${plot_ns}::w] <Shift-Down>     { ::cv_dashboard::zoom 16 }
    bind [set ${plot_ns}::w] <v>              { ::cv_dashboard::fit_vertically }
    bind [set ${plot_ns}::w] <h>              { ::cv_dashboard::fit_horizontally }
  }

  # Update frame to display frame marker in new plot
  update_frame internal $::cv_dashboard::mol w
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

  # Round to nearest frame number
  animate goto [expr { round(($x - $xplotmin) / $scalex + $xmin)}]
  if { $::cv_dashboard::track_frame == 0 } {
    # frame change doesn't trigger refresh, so we refresh manually
    refresh_values
  }
}


# Callback for click on marker in 2 cv plot
proc ::cv_dashboard::marker_clicked { index x y color marker } {

  animate goto [expr {$index - 1 }]
  if { $::cv_dashboard::track_frame == 0 } {
    # frame change doesn't trigger refresh, so we refresh manually
    refresh_values
  }
}


# Change frame in reaction to user input (arrow keys)
proc ::cv_dashboard::chg_frame { shift } {

  set nf [molinfo $::cv_dashboard::mol get numframes]

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
    refresh_values
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
  set max_f [expr [molinfo $::cv_dashboard::mol get numframes] - 1]
  if {$fmax > $max_f} { set fmax $max_f }

  $plothandle configure -xmin $fmin -xmax $fmax -plot
  update_frame internal $::cv_dashboard::mol w
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
  update_frame internal $::cv_dashboard::mol w
}


# Fit all x values within horizontal range, then fit vertically
proc ::cv_dashboard::fit_horizontally {} {
  variable ::cv_dashboard::plothandle

  $plothandle configure -xmin auto -xmax auto -ymin auto -ymax auto -plot
  update_frame internal $::cv_dashboard::mol w
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
