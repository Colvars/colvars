# Tcl package index file

# This file is generated *manually*
# and sourced either when an application starts up or
# by a "package unknown" script.  It invokes the
# "package ifneeded" command to set up package-related
# information so that packages will be loaded automatically
# in response to "package require" commands.  When this
# script is sourced, the variable $dir must contain the
# full path name of this file's directory.

proc init {} {
  set dir [file dirname [info script]]
  source [file join $dir cv_dashboard_update.tcl]
  set VERSION [::cv_dashboard::read_version $dir]
  set ::cv_dashboard::version $VERSION

  # Compare with version with local package
  set local_dir [file join [::cv_dashboard::get_local_dir] "cv_dashboard"]

  if { $local_dir ne [file dirname $dir] && [file exists $local_dir]} {
    # There is a local package separate from this instance
    # Compare versions
    set LOCAL_VERSION [::cv_dashboard::read_version $local_dir]

    if { [string compare $LOCAL_VERSION $VERSION] > 0 } {
      puts "Overriding current version ($VERSION) with local installation of cv_dashboard from $local_dir ($LOCAL_VERSION)"
      source [file join $local_dir pkgIndex.tcl]
      return
    }
  }

  # Now try to update local package from remote repository
  if [::cv_dashboard::self_update] {
    # Reload if necessary
    source [file join $local_dir pkgIndex.tcl]
    return
  }

  # Finish initializing this version
  package ifneeded cv_dashboard $VERSION "set env(CV_DASHBOARD_DIR) [list $dir]; [list source [file join $dir cv_dashboard.tcl]]"
}

init
