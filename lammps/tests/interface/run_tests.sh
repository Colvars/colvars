#!/bin/sh
if [ $# -lt 1 ] || [ ! -x "$1" ]
then
    cat <<EOF
 usage: $0 <path-to-lammps>

EOF
exit 1
fi
lmp="$@"
inputs=$(echo in.[0-9][0-9]? | sed -e s/in\.//g)

for id in $inputs
do \
    log="log.${id}"
    ${lmp} -echo screen -log none -in in.${id} \
      | egrep -v '(LAMMPS|OpenMP|MPI|serial|Memory|Loop|Perform|FFT|version)' \
      | sed -e 's/-0\.0000000000/0.0000000000 /g' -e '/^Section.*/,$d' \
            -e 's,\.cpp:[0-9]\+),.cpp),' > ${log}
    test -f ref/${log} || touch ref/${log}
    cmp ref/${log} ${log} && rm ${log} || diff -u ref/${log} ${log}
    rm -f out.colvars.state out.colvars.traj
    for f in ${id}.colvars.traj ${id}.colvars.state
    do \
        if [ -f $f ]
        then \
            test -f ref/$f || touch ref/$f
            cmp ref/$f $f && rm $f || diff -u ref/$f $f
        fi
    done
done

../library/run_tests.sh "$@"
