if declare -f module >& /dev/null ; then
    # Default modulefile in OpenMPI RPM package
    if module load mpi/openmpi-x86_64 >& /dev/null ; then
        echo "Loaded OpenMPI version $(mpirun --version | head -n 1 | rev | cut -d' ' -f1 | rev)"
    fi
fi
