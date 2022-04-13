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
  if { $type == "histogram" && $total_dim != 1 } {
    tk_messageBox -icon error -title "Colvars Dashboard Error"\
      -message "Select exactly 1 scalar quantity for a histogram plot.\n"
    return
  }

  set nf [molinfo $::cv_dashboard::mol get numframes]
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
  } elseif { $type == "histogram"} {
    set xname [lindex $name_list 0]
    # Save list of values for navigating
    set ::cv_dashboard::histogram_time_series $y($xname)
    set ::cv_dashboard::histogram_sorted_frames [lsort -real -indices $y($xname)]
    set ::cv_dashboard::histogram_sorted_values [list]

    set i 0
    foreach f $::cv_dashboard::histogram_sorted_frames {
      set ::cv_dashboard::histogram_frame_rank($f) $i
      lappend ::cv_dashboard::histogram_sorted_values [lindex $::cv_dashboard::histogram_time_series $f]
      incr i
    }

    lassign [compute_histogram $::cv_dashboard::histogram_sorted_values] delta centers frequencies

    set nbins [llength $centers]
    if { $nbins == 0 } { return }

    # Create plot with real freq data to enable exports
    # do not display lines but call plot to compute sizes
    set plothandle [multiplot -title "Histogram for colvar $cvs \[click, keyb arrows (+ Shift/Ctrl) to navigate\]" \
      -xlabel $xname -ylabel "N" -nostats \
      -xmin [expr [lindex $centers 0] - (0.5*$delta)] -xmax [expr [lindex $centers end] + (0.5*$delta)] \
      -ymin 0.0 -x $centers -y $frequencies -nolines -plot]

    set ns [namespace qualifiers $plothandle]
    # force bars to start at zero
    set ymin 0.0

    for {set j 0} {$j < $nbins} {incr j} {
      set left [expr [lindex $centers $j] - (0.5 * $delta)]
      set right [expr [lindex $centers $j] + (0.5 * $delta)]
      $plothandle draw rectangle $left $ymin $right [lindex $frequencies $j] -fill "#ecf6ff" -tags rect$j
    }

    set maxfreq 0
    foreach f $frequencies {
      if { $f > $maxfreq } { set maxfreq $f }
    }

    # Compute and plot cumulative distribution
    set delta [expr {$maxfreq * 1.0 / $nf}] ;# 1.0 to force floating-point operation
    set cumul_dist_x [list]
    set cumul_dist_y [list]
    set npoints 200 ;# target number of points (plot will have at most twice that number)
    set stride [expr {$nf < $npoints ? 1 : int($nf / $npoints)}]
    for { set i 0 } { $i < $nf } { incr i $stride } {
      lappend cumul_dist_x [lindex $::cv_dashboard::histogram_sorted_values $i]
      lappend cumul_dist_y [expr {($i + 1) * $delta}]
    }
    # Add missing final point if stride is not a divisor of nf
    if {[lindex $cumul_dist_y end] < $maxfreq } {
      lappend cumul_dist_x [lindex $::cv_dashboard::histogram_sorted_values end]
      lappend cumul_dist_y $maxfreq
    }
    $plothandle add $cumul_dist_x $cumul_dist_y -linecolor black -legend "Cumulative distribution"
  }

  $plothandle replot
  # bind mouse and keyboard events to callbacks
  set plot_ns [$plothandle namespace]
  set w [set ${plot_ns}::w]

  if { $type == "timeline" } {
    traj_animation_bindings $w
    bind $w <Button-1>       { ::cv_dashboard::plot_clicked %x %y }
    bind $w <Up>             { ::cv_dashboard::zoom 0.25 }
    bind $w <Down>           { ::cv_dashboard::zoom 4 }
    bind $w <Shift-Up>       { ::cv_dashboard::zoom 0.0625 }
    bind $w <Shift-Down>     { ::cv_dashboard::zoom 16 }
    bind $w <v>              { ::cv_dashboard::fit_vertically }
    bind $w <h>              { ::cv_dashboard::fit_horizontally }
  } elseif { $type == "2cv" } {
    traj_animation_bindings $w
  } elseif { $type == "histogram" } {
    bind $w <Button-1>       { ::cv_dashboard::plot_clicked %x %y }
    bind $w <Left>           { ::cv_dashboard::chg_frame -1 sorted }
    bind $w <Right>          { ::cv_dashboard::chg_frame 1 sorted }
    bind $w <Shift-Left>     { ::cv_dashboard::chg_frame -10 sorted }
    bind $w <Shift-Right>    { ::cv_dashboard::chg_frame 10 sorted }
    bind $w <Control-Left>   { ::cv_dashboard::chg_frame -50 sorted }
    bind $w <Control-Right>  { ::cv_dashboard::chg_frame 50 sorted }
    bind $w <Home>           { ::cv_dashboard::chg_frame start sorted }
    bind $w <End>            { ::cv_dashboard::chg_frame end sorted }
  }

  # Update frame to display frame marker in new plot
  update_frame internal $::cv_dashboard::mol w
}

# Create plot window for energy of biases
proc ::cv_dashboard::plot_bias_energy { } {
  variable ::cv_dashboard::plothandle
  set ::cv_dashboard::plottype timeline

  # Remove existing plot, if any
  if { [info exists plothandle] } {
    catch {$plothandle quit}
    unset plothandle
  }

  set biases [selected_biases]
  if { [llength $biases] == 0 } {
    # If no selection, plot all variables
    set biases [run_cv list biases]
    if { [llength $biases] == 0 } {
      return
    }
  }

  foreach b $biases {
    set y($b) [list]
  }

  set nf [molinfo $::cv_dashboard::mol get numframes]
  # Get list of values for all frames
  for {set f 0} {$f< $nf} {incr f} {
    run_cv frame $f
    run_cv update
    foreach b $biases {
      set val [run_cv bias $b update]
      # Catch NaN energies
      if { $val != $val } { set val 0.0 }
      lappend y($b) $val
    }
  }

  set plothandle [multiplot \
    -title {Bias energy trajectory   [click, keyb arrows (+ Shift/Ctrl) to navigate & zoom, v/h to fit vert/horizontally]} \
    -xlabel "Frame" -ylabel "Value" -nostats]
  set x {}
  for {set f 0} {$f < $nf} {incr f} { lappend x $f }
  foreach n $biases {
    $plothandle add $x $y($n) -legend $n
  }

  $plothandle replot
  # bind mouse and keyboard events to callbacks
  set plot_ns [$plothandle namespace]
  set w [set ${plot_ns}::w]

  traj_animation_bindings $w
  bind $w <Button-1>       { ::cv_dashboard::plot_clicked %x %y }
  bind $w <Up>             { ::cv_dashboard::zoom 0.25 }
  bind $w <Down>           { ::cv_dashboard::zoom 4 }
  bind $w <Shift-Up>       { ::cv_dashboard::zoom 0.0625 }
  bind $w <Shift-Down>     { ::cv_dashboard::zoom 16 }
  bind $w <v>              { ::cv_dashboard::fit_vertically }
  bind $w <h>              { ::cv_dashboard::fit_horizontally }


  # Update frame to display frame marker in new plot
  update_frame internal $::cv_dashboard::mol w
}


# Callback for click inside plot window, at coords x y
proc ::cv_dashboard::plot_clicked { x y } {

  set ns [$::cv_dashboard::plothandle namespace]
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

  if { $::cv_dashboard::plottype == "timeline" } {
    # Round to nearest frame number
    set newframe [expr {round(($x - $xplotmin) / $scalex + $xmin)}]
  } elseif { $::cv_dashboard::plottype == "histogram" } {
    set val [expr {($x - $xplotmin) / $scalex + $xmin}]
    # Find frame with closest value
    set newframe [lindex $::cv_dashboard::histogram_sorted_frames [bisect $::cv_dashboard::histogram_sorted_values $val]]
  }
  animate goto $newframe
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
proc ::cv_dashboard::chg_frame { shift { order "traj" } } {

  set nf [molinfo $::cv_dashboard::mol get numframes]

  if { $shift == "start" } {
    set f 0
  } elseif { $shift == "end" } {
    set f [expr $nf - 1]
  } else {
    set f [expr $::cv_dashboard::current_frame]
    if { $order == "sorted" } {
      # Work with frame rank instead of frame ID
      set f [expr $::cv_dashboard::histogram_frame_rank($f) + $shift]
    } else {
      set f [expr $f + $shift]
    }
  }

  # Keep within bounds [[O, numframes[[
  if { $f < 0 } { set f 0 }
  if { $f >= $nf } { set f [expr $nf - 1] }
  if { $order == "sorted" } {
    # Convert back to frame ID
    set f [lindex $::cv_dashboard::histogram_sorted_frames $f]
  }

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
  set timeline_color "blue"
  set histogram_color "red"
  set 2cv_color "red"

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
        $canv create line  $x $y1 $x $y2 -fill $timeline_color -tags frame_marker
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
        $canv create oval $x1 $y1 $x2 $y2 -outline white -fill $2cv_color -tags frame_marker
      } elseif { $::cv_dashboard::plottype == "histogram" } {
        set xmin [set ${ns}::xmin]
        set xmax [set ${ns}::xmax]

        set y1 [set ${ns}::yplotmin]
        set y2 [set ${ns}::yplotmax]
        set v [lindex $cv_dashboard::histogram_time_series $f]

        set xplotmin [set ${ns}::xplotmin]
        set scalex [set ${ns}::scalex]
        set x [expr $xplotmin+($scalex*($v-$xmin))]

        set canv "[set ${ns}::w].f.cf"
        $canv delete frame_marker
        $canv create line  $x $y1 $x $y2 -fill $histogram_color -tags frame_marker
      }
    }
  }
}


# Create plot window for energy of biases
proc ::cv_dashboard::compute_histogram { sorted_values } {
  set nbins $::cv_dashboard::nbins

  if {[llength $sorted_values] < 2} {
    tk_messageBox -icon error -title "Colvars Dashboard Error"\
      -message "At least two frames are necessary to compute a histogram.\n"
    return "" ""
  }
  set min [lindex $sorted_values 0]
  set max [lindex $sorted_values end]
  set delta [expr ($max - $min) / ($nbins - 1)]
  if { $delta == 0. } {
    set delta 1.
  } else {
    # Adjust to something round in decimal terms
    set delta [simplify $delta]
  }
  # Align bin boundaries on integer multiples of delta
  set min [expr floor($min / $delta) * $delta]
  # Adjust bins as required
  set nbins [expr int(($max - $min) / $delta) + 1]

  for {set i 0} {$i < $nbins} {incr i} {
    set c($i) 0
  }
  foreach v $sorted_values {
    incr c([expr {int(floor(($v-$min)/$delta))}])
  }
  set centers [list]
  set freqs [list]
  for {set i 0} {$i < $nbins} {incr i} {
    lappend centers [expr {$min + $delta * ($i + 0.5)}]
    lappend freqs $c($i)
  }
  return [list $delta $centers $freqs]
}


proc ::cv_dashboard::simplify x {
  #Compute nearest "round" number, with first 2 decimals:
  # 1, 1.5, 2, 2.5, 3, 4, 5, ..., 9

  # Factor has two components:
  set pow10 [expr 10**(floor(log10($x)))]  ;# -> bring x between 1 and 10
  set xn [expr $x / $pow10]
  set half_int [expr floor(log10(3.4 * $xn) + 1) / 2]  ;#-> use half-integer steps until .3, integer steps above

  set factor [expr $pow10 * $half_int]
  set x [expr round($x / $factor) * $factor]
  return $x
}
