export mol_name=ammonia_water
export numsteps=50000000
export run=$(basename ${PWD})

# Create folders
for dir in replica-{000..003} ; do
    mkdir -p $dir
    unset dir
done

if [ -z "${NAMD_HOME}" ] ; then
    NAMD_HOME=../NAMD_colvars-bias-exchange/Linux-x86_64-g++.netlrts
fi

for ijob in {0..1} ; do
    export ijob
    $NAMD_HOME/charmrun ++local +p4 $NAMD_HOME/namd2 \
        +replicas 4 ${run}.namd +stdout replica-%03d/${mol_name}.${run}.replica-%03d.$(printf "%04d" ${ijob}).out
    unset ijob
done

unset mol_name numsteps run
