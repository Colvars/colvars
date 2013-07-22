#!/bin/sh
if [ $# -lt 1 ] 
then
    cat <<EOF
 usage: $0 <command-to-run-lammps>

EOF
exit 1
fi

./clean.sh

log=log.screen
$@ -echo screen -log none -in in.peptide \
    | egrep -v '(LAMMPS|OpenMP|MPI|serial|Memory|Loop|Perform|FFT|version)' \
    | sed -e 's/-0\.0000000000/0.0000000000 /g' -e 's,\.cpp:[0-9]\+),.cpp),' \
    | tee ${log}

sed -i  -e 's/^Loop.*$//' -e 's/^Memory.*$//' \
	-e 's/-0\.0000000000/0.0000000000 /g' \
	-e 's/\.cpp:[0-9]\+)/.cpp)/' log.*-?

# compare fix spring to fix colvars
# these always show some small differences.
for s in 1 2 3 4 5
do \
    sed -i -e '/^colvars/d' log.colvars-$s
    diff -u log.spring-$s log.colvars-$s
done

# do an inexact compare on colvars trajectory files
deleteme=
for r in $(/bin/ls reference/ | grep .traj)
do \
    test -f reference/$r || touch reference/$r
    ../common/cvtraj-compare.pl reference/$r $r && deleteme="${deleteme} $r"
done

rm -f $deleteme
