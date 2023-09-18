# Tcl package index file, version 1.1

# This file is generated manually
# and sourced either when an application starts up or
# by a "package unknown" script.  It invokes the
# "package ifneeded" command to set up package-related
# information so that packages will be loaded automatically
# in response to "package require" commands.  When this
# script is sourced, the variable $dir must contain the
# full path name of this file's directory.

set dir [file dirname [info script]]
set version_file [open "${dir}/VERSION"]
gets $version_file CV_DASHBOARD_VERSION
close $version_file
# Convert to Tcl-style package version number
set CV_DASHBOARD_TCL_VERSION [string map { "-" "." } $CV_DASHBOARD_VERSION]

package ifneeded cv_dashboard $CV_DASHBOARD_TCL_VERSION "set env(CV_DASHBOARD_DIR) [list $dir]; [list source [file join $dir cv_dashboard.tcl]]"
