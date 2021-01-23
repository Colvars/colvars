# Performs Adiabatic Biased MD in Colvars / NAMD
# after Marchi & Ballone JCP 1999 / Paci & Karplus JMB 1999
# (not to be confused with Adaptively Biased MD)
#
# Jérôme Hénin <henin@ibpc.fr>

# Usage:
#
# source abmd.tcl
# setup_ABMD x 10. 42.  ;# colvar k x_max
#
# cv config "
#   scriptedColvarForces on
#   colvar {
#     name x
# ..."

proc setup_ABMD { colvar force_k cv_max } {
  namespace eval ::ABMD {}
  set ::ABMD::cvname $colvar
  set ::ABMD::xmax $cv_max
  set ::ABMD::k $force_k
}

proc calc_colvar_forces { ts } {
  if { $ts == 0 } {
    set ::ABMD::max [cv colvar $::ABMD::cvname value]
  }
  namespace eval ::ABMD {
    set x [cv colvar $cvname value]
    if { $x > $max } {
      if { $x <= $xmax } { set max $x }
    } else {
      cv colvar $cvname addforce [expr { $k * ($max - $x) } ]
    }
  }
}

