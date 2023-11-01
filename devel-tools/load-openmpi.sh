# Ensure that Lmod is properly initialized
source /etc/profile

if declare -f module >& /dev/null ; then
    # Default modulefile in OpenMPI RPM package
    if module load mpi/openmpi-x86_64 >& /dev/null ; then
        echo "Loaded OpenMPI version $(mpirun --version | head -n 1 | rev | cut -d' ' -f1 | rev)"
        if ompi_info --param btl vader | grep -q vader ; then
            # Enforce intra-node communication on recent OpenMPI versions
            export OMPI_MCA_btl="vader,self"
        fi
    fi
fi
