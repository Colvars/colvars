####################################
# Self-update utilities

namespace eval ::cv_dashboard {
  variable repo_user "Colvars"
  variable repo_name "colvars"
  variable branch "dashboard_update"
}

proc ::cv_dashboard::check_version_file {} {
  variable version
  variable repo_user
  variable repo_name
  variable branch

  # Download VERSION file from repo
  set version_url "https://raw.githubusercontent.com/$repo_user/$repo_name/$branch/vmd/cv_dashboard/VERSION"
  set temp_dir [get_local_dir]
  file mkdir $temp_dir
  set temp_file [file join $temp_dir "remote_version.txt"]

  if {[catch {exec curl -L -s -o $temp_file $version_url}] &&
      [catch {exec wget -O $temp_file $version_url}]} {
    puts "exec curl -L -s -o $temp_file $version_url"
    puts "exec wget -q -O $temp_file $version_url"
    error "Could not check version file"
  }

  set fp [open $temp_file r]
  set remote_version [string map { "-" "." } [string trim [read $fp]]]
  close $fp
  file delete $temp_file

  set comparison [compare_versions $version $remote_version]
  return [list $comparison $remote_version]
}

proc ::cv_dashboard::compare_versions {current latest} {
  # Versioned by date: e.g. 2025-05-21
  set current_parts [split $current "."]
  set latest_parts [split $latest "."]

  # Pad with zeros if needed
  while {[llength $current_parts] < 3} {lappend current_parts 0}
  while {[llength $latest_parts] < 3} {lappend latest_parts 0}

  for {set i 0} {$i < 3} {incr i} {
    set curr [lindex $current_parts $i]
    set lat [lindex $latest_parts $i]

    if {$lat > $curr} {
      return "outdated"
    } elseif {$lat < $curr} {
      return "newer"
    }
  }
  return "current"
}

proc ::cv_dashboard::download_github_directory {user repo directory {branch "main"} {destdir "."}} {
  set url "https://github.com/$user/$repo/archive/refs/heads/$branch.zip"
  set archive "[file join $destdir $repo.zip]"
  set temp_dir "[file join $destdir temp_$repo]"

  file mkdir $destdir
  file mkdir $temp_dir

  puts "Downloading $url"
  # Download full repo as ZIP
  set downloaded 0
  if {![catch {exec curl -L -s -o $archive $url}]} {
    set downloaded 1
  } elseif {![catch {exec wget -q -O $archive $url}]} {
    set downloaded 1
  } elseif {$::tcl_platform(platform) eq "windows" &&
          ![catch {exec powershell -command "Invoke-WebRequest -Uri '$url' -OutFile '$archive'" 2>nul}]} {
    set downloaded 1
  }

  if {!$downloaded} {
    error "Failed to download the repository archive."
  }

  # Extract to temp directory
  set extracted 0

  puts "Extracting archive..."
  if {![catch {exec unzip -q $archive -d $temp_dir}]} {
    set extracted 1
  } elseif {$::tcl_platform(platform) eq "windows"} {
    set ps_extract_script {
      Add-Type -AssemblyName System.IO.Compression.FileSystem
      [System.IO.Compression.ZipFile]::ExtractToDirectory('%ARCHIVE%', '%DEST%')
    }
    set ps_script [string map [list %ARCHIVE% $archive %DEST% $temp_dir] $ps_extract_script]
    if {![catch {exec powershell -NoProfile -Command $ps_script}]} {
      set extracted 1
    }
  }

  if {!$extracted} {
    error "Failed to extract ZIP archive."
  }
  file delete $archive

  # Move specific directory to destination
  set extracted_repo "[file join $temp_dir $repo-$branch]"
  set source_dir "[file join $extracted_repo $directory]"
  set dest_dir "[file join $destdir [file tail $directory]]"

  if {![file exists $source_dir]} {
    file delete -force $temp_dir
    error "Directory '$directory' not found in repository"
  }

  file rename $source_dir $dest_dir
  file delete -force $temp_dir

  puts "Directory downloaded to: $dest_dir"
  return $dest_dir
}

proc ::cv_dashboard::download_github_directory_git {user repo directory {branch "main"} {destdir "."}} {
  set repo_url "https://github.com/$user/$repo.git"
  set clone_dir "[file join $destdir $repo-sparse]"

  file mkdir $destdir
  # Check if git is available
  if {[catch {exec git --version 2>/dev/null}]} {
    error "Git not available - use download_github_directory instead"
  }

  # Clone with no checkout
  exec git clone --no-checkout --depth 1 --branch $branch $repo_url $clone_dir 2>@1

  # Set up sparse checkout
  set git_dir [file join $clone_dir ".git"]
  exec git --git-dir=$git_dir --work-tree=$clone_dir config core.sparseCheckout true

  # Specify directory to checkout
  set sparse_file [file join $clone_dir ".git" "info" "sparse-checkout"]
  set fp [open $sparse_file w]
  puts $fp "$directory/*"
  close $fp

  # Checkout only the specified directory
  exec git --git-dir=$git_dir --work-tree=$clone_dir checkout 2>@1

  # Move the directory to final location
  set source_dir "[file join $clone_dir $directory]"
  set dest_dir "[file join $destdir [file tail $directory]]"

  if {![file exists $source_dir]} {
    file delete -force $clone_dir
    error "Directory '$directory' not found in repository"
  }

  file rename $source_dir $dest_dir
  file delete -force $clone_dir

  puts "Directory downloaded to: $dest_dir"
  return $dest_dir
}

proc ::cv_dashboard::download_directory_smart {user repo directory {branch "main"} {destdir "."}} {
  # Try git sparse-checkout first (most efficient)
  if {![catch {exec git --version 2>/dev/null}]} {
    puts "Using git sparse-checkout..."
    return [download_github_directory_git $user $repo $directory $branch $destdir]
  } else {
    puts "Git not available, downloading full repository..."
    return [download_github_directory $user $repo $directory $branch $destdir]
  }
}

proc ::cv_dashboard::self_update {{force false}} {
  variable version
  variable repo_user
  variable repo_name
  variable branch

  # Check for updates
  if {!$force} {
    lassign [::cv_dashboard::check_version_file] status latest_version
    if {$status eq "current"} {
      puts "Local cv_dashboard ($version) is up to date"
      return false
    } elseif {$status eq "newer"} {
      puts "Local cv_dashboard ($version) is newer than remote ($latest_version)"
      return false
    }
    puts "cv_dashboard update available: $version -> $latest_version"
  }

  set package_dir [get_local_dir]
  set backup_dir "${package_dir}.backup.[clock seconds]"

  # Create backup
  if {[file exists $package_dir]} {
    file copy $package_dir $backup_dir
  }

  try {
    # Download new version
    puts "Downloading update..."
    set temp_dir [file join [file dirname $package_dir] "temp_update"]
    download_directory_smart $repo_user $repo_name "vmd/cv_dashboard" $branch $temp_dir

    # Replace current installation
    file delete -force $package_dir
    file rename $temp_dir $package_dir

    puts "Update completed successfully"
    if {[file exists $backup_dir]} {
      puts "Backup available at: $backup_dir"
    }
    cleanup_backups $package_dir
    return true

  } on error {err} {
    puts "Update failed: $err"

    if {[file exists $backup_dir]} {
      puts "Restoring from backup..."
      if {[file exists $package_dir]} {
        file delete -force $package_dir
      }
      file rename $backup_dir $package_dir
    }
    return false
  }
}

proc ::cv_dashboard::get_local_dir {} {
  set dirname "colvars"
  switch $::tcl_platform(platform) {
    "windows" {
      set base [expr {[info exists ::env(LOCALAPPDATA)] ? $::env(LOCALAPPDATA) : $::env(USERPROFILE)}]
      return [file join $base $dirname]
    }
    "unix" {
      set home $::env(HOME)
      if {$::tcl_platform(os) eq "Darwin"} {
        return [file join $home "Library" "Application Support" $dirname]
      } else {
        return [file join $home ".local" "share" $dirname]
      }
    }
  }
}

proc ::cv_dashboard::cleanup_backups { package_dir {keep 2}} {
  set parent_dir [file dirname $package_dir]
  set package_name [file tail $package_dir]

  # Find backup directories
  set backups {}
  foreach item [glob -nocomplain -directory $parent_dir "${package_name}.backup.*"] {
    if {[file isdirectory $item]} {
      lappend backups $item
    }
  }

  # Sort by timestamp and remove old ones
  set backups [lsort $backups]
  set to_remove [lrange $backups 0 end-$keep]

  foreach backup $to_remove {
    puts "Removing old backup: $backup"
    file delete -force $backup
  }
}

proc ::cv_dashboard::read_version { path } {
  set version_file [open [file join $path VERSION]]
  gets $version_file git_date
  set version [string map { "-" "." } $git_date]
  close $version_file
  return $version
}
