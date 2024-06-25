#!/bin/sh

# Usage: ./calc_colvars.sh config.colvars structure.<psf/gro> traj.<dcd/xtc/trr> [...]

echo "Writing trajectory for configuration $1 to colvars.traj"

echo "cv molid 0; cv configfile $1" > _calc_cvs.tmp.tcl

cat >> _calc_cvs.tmp.tcl << EOF
set fileName "colvars.traj"
set o [open \$fileName w]
puts -nonewline \$o [cv printframelabels]
set nf [molinfo top get numframes]
for {set f 0} {\$f< \$nf} {incr f} {
  cv frame \$f
  cv update
  puts -nonewline \$o [cv printframe]
}
close \$o
exit
EOF

shift
vmd -dispdev text $* -e _calc_cvs.tmp.tcl > vmd.log
rm -f _calc_cvs.tmp.tcl

