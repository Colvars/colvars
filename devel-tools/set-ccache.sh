if [ -d /usr/lib/ccache ] ; then
    export CCACHE_HOME=/usr/lib/ccache
fi
if [ -d /usr/lib64/ccache ] ; then
    export CCACHE_HOME=/usr/lib64/ccache
fi

export PATH=${CCACHE_HOME}:${PATH}

if [ -z "${CCACHE_DIR}" ] || [ "x${CCACHE_DIR}" == "x/var/cache/ccache" ] ; then
    # Test for default setting in RHEL-like distributions
    if grep -sq 'CCACHE_DIR=/var/cache/ccache' /etc/profile.d/ccache.sh ; then
        if ! mktemp /var/cache/ccache/test.XXXXXX >& /dev/null ; then
            # Prevent using unwritable path inside containers
            export CCACHE_DIR=~/.ccache
        fi
    fi
fi

# Report cache statistics
ccache -s
