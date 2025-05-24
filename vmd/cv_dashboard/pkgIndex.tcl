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
  set version_file [open "${dir}/VERSION"]
  gets $version_file git_date
  close $version_file
  # Convert to Tcl-style package version number
  set CV_DASHBOARD_VERSION [string map { "-" "." } $git_date]

  # Compare with version with local package
  source [file join $dir cv_dashboard_update.tcl]
  set package_dir [::cv_dashboard::get_local_dir]

  if { $package_dir ne $dir && [file exists $package_dir]} {
    # Compare versions
    set version_file [open [file join $package_dir "cv_dashboard" VERSION]]
    gets $version_file git_date
    close $version_file
    set LOCAL_VERSION [string map { "-" "." } $git_date]

    if { [string compare $LOCAL_VERSION $CV_DASHBOARD_VERSION] > 0 } {
      puts "Overriding current version ($CV_DASHBOARD_VERSION) with local installation of cv_dashboard from $package_dir ($LOCAL_VERSION)"
      source [file join $package_dir "cv_dashboard" pkgIndex.tcl]
      return
    }
  }

  # Only update local install
  if { $package_dir eq $dir } {
    set updated [::cv_dashboard::self_update]

    if $updated {
      # Source ourselves! Hoping this doesn't loop
      source [file join $package_dir "cv_dashboard" pkgIndex.tcl]
      return
    }
  }

  # Finish initializing package
  package ifneeded cv_dashboard $CV_DASHBOARD_VERSION "set env(CV_DASHBOARD_DIR) [list $dir]; [list source [file join $dir cv_dashboard.tcl]]"
}

init
