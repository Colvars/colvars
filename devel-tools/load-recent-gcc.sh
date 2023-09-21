# Load Software Collections package for RHEL/CentOS 7
if [ -s /etc/redhat-release ] && [ -s /opt/rh/devtoolset-10/enable ] ; then
    source /opt/rh/devtoolset-10/enable && echo "Loaded: $(gcc --version)" | head -n 1
fi
