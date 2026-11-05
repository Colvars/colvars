FROM quay.io/centos/centos:stream9

LABEL org.opencontainers.image.description="Development environment for CentOS Stream Linux 9 with Intel oneAPI installed"

SHELL ["/bin/bash", "-c"]

RUN dnf -y update && \
    dnf -y install dnf-plugins-core && \
    dnf config-manager --set-enabled crb && \
    dnf -y install \
        https://dl.fedoraproject.org/pub/epel/epel-release-latest-9.noarch.rpm \
        https://dl.fedoraproject.org/pub/epel/epel-next-release-latest-9.noarch.rpm && \
    dnf -y install \
        "@Development Tools" \
        gcc gcc-c++ \
        git cmake ccache ninja-build mold \
        ncurses-devel \
        Lmod openmpi-devel tbb-devel fftw-devel tcl-devel && \
    dnf clean all

# Alternatively, install the full intel-basekit and intel-hpckit
RUN printf '%s\n' \
        '[oneAPI]' \
        'name=Intel® oneAPI repository' \
        'baseurl=https://yum.repos.intel.com/oneapi' \
        'enabled=1' \
        'gpgcheck=1' \
        'repo_gpgcheck=1' \
        'gpgkey=https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB' \
        > /etc/yum.repos.d/oneAPI.repo && \
    dnf -y install intel-oneapi-compiler-dpcpp-cpp intel-oneapi-mpi-devel && \
    dnf clean all

RUN printf '%s\n' \
        'if [ -z "${ONEAPI_ROOT}" ] ; then' \
        '    echo "Usage of Intel oneAPI is subject to the license agreement at: https://www.intel.com/content/www/us/en/developer/articles/license/end-user-license-agreement.html" >&2' \
        '    source /opt/intel/oneapi/setvars.sh > /dev/null' \
        'fi' > /etc/profile.d/oneapi.sh && \
    chmod 644 /etc/profile.d/oneapi.sh

CMD ["/bin/bash", "-l"]
