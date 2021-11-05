# Load Software Collections packages on a RedHat-style distribution
if hash scl >& /dev/null ; then
    if scl -l | grep -q ^devtoolset-10 ; then
        source /opt/rh/devtoolset-10/enable && \
            echo "Loaded: $(gcc --version)"
    fi
fi
