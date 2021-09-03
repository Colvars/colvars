if [ -d /usr/lib/ccache ] ; then
    export CCACHE_HOME=/usr/lib/ccache
fi
if [ -d /usr/lib64/ccache ] ; then
    export CCACHE_HOME=/usr/lib64/ccache
fi

export PATH=${CCACHE_HOME}:${PATH}

# Report cache statistics
ccache -s
