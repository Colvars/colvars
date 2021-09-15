# Load Software Collections packages on a RedHat-style distribution
if hash scl >& /dev/null ; then
    if scl -l | grep -q ^rh-git227 ; then
        source /opt/rh/rh-git227/enable && \
            echo "Loaded: $(git --version)"
    fi
fi
