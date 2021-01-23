# Performs Adiabatic Biased MD in Colvars / NAMD
# after Marchi & Ballone JCP 1999 / Paci & Karplus JMB 1999
# (not to be confused with Adaptively Biased MD)
#
# Jérôme Hénin <henin@ibpc.fr>

# Usage:
#
# source abmd.tcl
# setup_ABMD x 10. 42. [up/down] ;# colvar k z_stop
#
# cv config "
#   scriptedColvarForces on
#   colvar {
#     name x
# ..."

# Calculates the high-water mark z_max of colvar cvname
# Ignores its arguments (dummy component)
# Note: this is not restartable (z_max info is lost at the end of a run)

proc calc_z_max { args } {
  namespace eval ::ABMD {
    set z [cv colvar $cvname value]
    if { ![info exists z_max] } {
      set z_max $z
      puts "ABMD: z_max = $z_max"
    } elseif { ($z_max < $z_stop) && ($z >= $z_stop) } {
      set z_max $z_stop
      puts "ABMD: z_max = $z_max"
    } elseif { ($z > $z_max) && ($z <= $z_stop) } {
      set z_max $z
      puts "ABMD: z_max = $z_max"
    }
    return $z_max
  }
}

proc calc_z_min { args } {
  namespace eval ::ABMD {
    set z [cv colvar $cvname value]
    if { ![info exists z_min] } {
      set z_min $z
      puts "ABMD: z_min = $z_min"
    } elseif { ($z_min > $z_stop) && ($z <= $z_stop) } {
      set z_min $z_stop
      puts "ABMD: z_min = $z_min"
    } elseif { ($z < $z_min) && ($z >= $z_stop) } {
      set z_min $z
      puts "ABMD: z_min = $z_min"
    }
    return $z_min
  }
}

proc calc_z_max_gradient { args } { return 0 }
proc calc_z_min_gradient { args } { return 0 }

proc setup_ABMD { colvar force_k z_stop {direction up} } {
  # cv config "scriptedColvarForces on"

  namespace eval ::ABMD {}
  set ::ABMD::cvname $colvar
  set ::ABMD::z_stop $z_stop
  set ::ABMD::k $force_k
  set ::ABMD::direction $direction

# Variable with dummy component - all the work is done by scripted function
  if { $direction == "up" } {
    set functionName "z_max"
  } elseif { $direction == "down" } {
    set functionName "z_min"
  } else {
    puts "ABMD ERROR: could not parse direction \"$direction\" -- should be up or down"
    puts ""
    exit 1
  }
  cv config "
colvar {
name z_ref
  distance {
    group1 { dummyAtom (0, 0, 0) }
    group2 { dummyAtom (0, 0, 0) }
  }
  scriptedFunction $functionName
}
"
}


proc calc_colvar_forces { ts } {
  namespace eval ::ABMD {

    set z [cv colvar $cvname value]
    set z_ref [cv colvar z_ref value]

    # We need to test because we could be beyond z_stop and should not apply a force then
    # even if z != z_ref
    if { ($direction == "up" && $z < $z_ref) || ($direction == "down" && $z > $z_ref)} {
      cv colvar $cvname addforce [expr { $k * ($z_ref - $z) } ]
    }
  }
}

