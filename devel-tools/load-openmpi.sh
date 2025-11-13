# Ensure that the module command is properly initialized
if [ -f /etc/profile.d/modules.sh ] ; then
    source /etc/profile.d/modules.sh
fi

if declare -f module >& /dev/null && [ -d /usr/share/modulefiles ] ; then
    module use /usr/share/modulefiles
    # Default modulefile in OpenMPI RPM package
    if module load mpi/openmpi-x86_64 ; then
        echo "Loaded OpenMPI version $(mpirun --version | head -n 1 | rev | cut -d' ' -f1 | rev)"
        if ompi_info --param btl vader | grep -q vader ; then
            # Enforce intra-node communication on recent OpenMPI versions
            export OMPI_MCA_btl="vader,self"
        fi
    fi
fi
