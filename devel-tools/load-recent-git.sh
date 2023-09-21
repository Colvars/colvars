# Load Software Collections package for RHEL/CentOS 7
if [ -s /etc/redhat-release ] && [ -s /opt/rh/rh-git227/enable ] ; then
    source /opt/rh/rh-git227/enable && echo "Loaded: $(git --version)"
fi
