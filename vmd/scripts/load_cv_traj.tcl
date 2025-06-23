# Load cv trajectories
# To be used as scripted colvars
# Only one copy - no way to distinguish different molecules so far
# Ideally the data should be associated to a given VMD molecule

if { ![array exists ::colvar_trajectory ]} {
  array set ::colvar_trajectory {}
}


proc load_cv_traj { filenames {molid "top"}} {

  if { $molid == "top" } {
    set molid [molinfo top]
  }
  # Detect and parse glob patterns in input filenames
  if [regexp {\*|\?} $filenames] {
    set filenames [lsort [glob $filenames]]
  }

  set steps [list]
  set file_num 0
  set title ""

  # (1) read raw data as list of lines containing values

  foreach n $filenames {
    puts "Opening $n ..."
    set fd [open $n "r"]
    set line_in_file 0
    set step -1  ;# MD step read from first column
    incr file_num

    while { [gets $fd line] >= 0 } {
      if { [llength $line] == 0} continue
      if { [lindex $line 0] == "#" } {
        if { $title == "" } {
          set title $line
          continue
        } elseif { $title != $line } {
          puts "Error: read incompatible title lines:\n$title\n$line"
          return
        } else {
          # New title matches previous, nothing to do
          continue
        }
      }
      # This is a data line
      set prev_step $step
      set step [lindex $line 0]
      if { $step == $prev_step } { continue }
      incr line_in_file
      if { $line_in_file == 1 && $file_num > 1 } {
        # Skip redundant first data line in subsequent files
        continue
      }
      # Convert vector quantities to Tcl lists: ( 1.0 , 2.3 ) -> { 1.0 2.3 }
      set line [string map {"(" "{" ")" "}" "," " "} $line]
      lappend steps $line
    }
    close $fd
  }

  set nsteps [llength $steps]
  puts "Read $nsteps steps from colvars.traj files"


  # (2) Use heuristic to match number of samples to number of frames
  # assuming that the molecular trajectory frequency is a multiple of the colvars.traj frequency

  set nf [molinfo $molid get numframes]

  if { $nf > 0 && [expr {($nsteps - 1) % $nf}] == 0 } {
    # VMD trajectory is missing first step
    set stride [expr {($nsteps - 1) / $nf}]
    set skip_steps $stride
    puts "Assuming $stride colvars.traj steps per trajectory frame, skipping $stride initial steps."
  } elseif { $nf > 1 && [expr {($nsteps - 1) % ($nf - 1)}] == 0 } {
    set stride [expr {($nsteps - 1) / ($nf - 1)}]
    set skip_steps 0
    puts "Assuming $stride colvars.traj steps per trajectory frame."
  } else {
    puts "Error: cannot match $nsteps steps with $nf frames."
    return
  }

  # (3) Store trajectory data in persistent array

  if { ![info exists ::cv_dashboard::colvar_trajectory($molid)] } {
    set ::cv_dashboard::colvar_trajectory($molid) [dict create]
  }

  set cvs [lrange $title 2 end]
  foreach cv $cvs {
    # Overwrite any previous data for the same molid and cv
    dict set ::cv_dashboard::colvar_trajectory($molid) $cv [list]
  }

  for {set f 0} {$f < $nf} {incr f} {
    set line [lindex $steps [expr $skip_steps + $stride * $f]]
    if {[llength $line] != [expr [llength $cvs] + 1]} {
      puts "Wrong number of elements in line:\n$line\nfor cv list\n$cvs"
      return
    }
    for { set i 0 } { $i < [llength $cvs]} { incr i } {
      set cv [lindex $cvs $i]
      dict lappend ::cv_dashboard::colvar_trajectory($molid) $cv [lindex $line [expr $i + 1]]
    }
  }

  # (4) Create scripted variables

  foreach cv $cvs {
    create_traj_colvar $molid $cv
  }
}


#  Create a scripted colvar that returns values from a precomputed trajectory

proc create_traj_colvar { molid cv } {

  set traj [dict get $::cv_dashboard::colvar_trajectory($molid) $cv]
  set dim [llength [lindex $traj 0]]
  set null [lrepeat $dim 0]

  # Specialized proc for given cv name, which returns 0 if data or frame is missing
  proc calc_traj_$cv { x } "
    set f \[cv frame\]
    set molid \[cv molid\]
    if { !\[info exists ::cv_dashboard::colvar_trajectory(\$molid)\] ||
          !\[dict exists \$::cv_dashboard::colvar_trajectory(\$molid) $cv\]} {
      return $null
    }

    set traj \[dict get \$::cv_dashboard::colvar_trajectory(\$molid) $cv\]
    if { \$f < \[llength \$traj\] } {
      return \[lindex \$traj \$f\]
    } else {
      return $null
    }
  "

  catch "cv colvar traj_$cv delete"

  set configString "
  colvar {
    name traj_$cv
    scriptedFunction traj_$cv"

  if { $dim > 1 } {
    append configString "
    scriptedFunctionType vector
    scriptedFunctionVectorSize $dim"
  }

  append configString "
    distance {
      # The distance component is just a placeholder
      group1 {
        dummyatom (0, 0, 0)
      }
      group2 {
        dummyatom (0, 0, 0)
      }
    }
  }
  "
  cv config $configString
}


# Sort frames in current molecule associated with Colvars
# into frames binned by values of a cv

proc bin_frames_by_cv { cv xmin dx nbins } {

  set src_molid [cv molid]
  set src [atomselect $src_molid "all"]
  $src writepsf "_cv_tmp.psf"

  # Create one target molecule per bin
  for { set  bin 0 } { $bin < $nbins } { incr bin } {
    # Duplicate the topology using a temp file (maybe there is a more elegant way)
    set molid($bin) [mol new "_cv_tmp.psf"]
    # Molecule name: cv and bin center
    mol rename $molid($bin) "$cv [expr {$xmin + ($bin + 0.5) * $dx}]"
    # keep an atomselect object per target molecule
    set dest($bin) [atomselect $molid($bin) "all"]
  }
  mol top $src_molid


  set nf [molinfo $src_molid get numframes]
  for {set f 0} {$f < $nf} {incr f} {
    # Display progress for impatient users
    if {$nf > 20 && [expr {$f % ($nf/20)}] == 0} { puts stdout "frame $f"}

    cv frame $f
    set x [cv colvar $cv update]
    set bin [expr { int(floor(($x-$xmin)/$dx)) }]
    if { $bin < 0 || $bin >= $nbins } {
      # Outside of target range
      continue
    }

    # create frame
    animate dup $molid($bin)
    # copy coordinates
    $dest($bin) frame [molinfo $molid($bin) get numframes]
    $src frame $f
    $dest($bin) set {x y z} [$src get {x y z}]
  }
}