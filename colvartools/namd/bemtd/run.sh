export mol_name=ammonia_water
export numsteps=50000000

# git clone git@github.com:TMB-CSB/charmm-ff

# Create empty folders
for dir in replica-{000..003} ; do
    mkdir -p $dir
    rm -f ${dir}/*
    unset dir
done

NAMD_HOME=NAMD_colvars-bias-exchange/Linux-x86_64-g++.netlrts

for ijob in {0..1} ; do
    export ijob
    $NAMD_HOME/charmrun ++local +p4 $NAMD_HOME/namd2 \
        +replicas 4 bemtd.namd +stdout replica-%03d/${mol_name}.bemtd.replica-%03d.$(printf "%04d" ${ijob}).out
    unset ijob
done

unset numsteps mol_name NAMD_HOME

