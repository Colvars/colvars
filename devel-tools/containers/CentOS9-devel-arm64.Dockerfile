FROM --platform=linux/arm64 quay.io/centos/centos:stream9

LABEL org.opencontainers.image.description="Development environment for CentOS Stream Linux 9"

# Bash is needed for brace expansion, "source" and Lmod's "module" command
SHELL ["/bin/bash", "-c"]

RUN dnf -y update && \
    dnf -y install dnf-plugins-core && \
    dnf config-manager --set-enabled crb && \
    dnf -y install \
        https://dl.fedoraproject.org/pub/epel/epel-release-latest-9.noarch.rpm \
        https://dl.fedoraproject.org/pub/epel/epel-next-release-latest-9.noarch.rpm && \
    dnf -y install \
        bc man-db vim emacs tcsh which bash-completion rsync screen tmux \
        "@Development Tools" \
        gcc gcc-c++ gcc-gfortran cppcheck \
        gcc-toolset-{12,13} clang clang-tools-extra \
        autoconf automake cvs git subversion mercurial \
        cmake ccache ninja-build mold \
        Lmod openmpi-devel tbb-devel fftw-devel tcl-devel \
        python3-{devel,tkinter,virtualenv,numpy,scipy} \
        ncurses-devel \
        rlwrap \
        tk-devel fltk-devel \
        mesa-dri-drivers mesa-libGL-devel libglvnd-devel \
        libXi-devel libXinerama-devel libpng-devel \
        netcdf-devel hdf5-devel sqlite-devel expat-devel \
        doxygen \
        texlive-latex texlive-epstopdf texlive-pdftex texlive-tex4ht \
        texlive-graphics texlive-metafont texlive-collection-latexrecommended \
        texlive-collection-fontsrecommended \
        latexmk latex2html \
        inkscape gimp wget \
        squashfs-tools && \
    dnf clean all

# Build Charm++
RUN umask 022 && \
    source /etc/profile && \
    module load mpi && \
    git clone --single-branch --depth=1 -b v8.0.1 https://github.com/UIUC-PPL/charm.git /opt/charm && \
    cd /opt/charm && \
    export CCACHE_DIR=/tmp && \
    ./build charm++ mpi-linux-arm8 -j16 --with-production && \
    ./build charm++ multicore-linux-arm8 -j16 --with-production && \
    ./build charm++ netlrts-linux-arm8 -j16 --with-production

# Build and install spiff
RUN git clone https://github.com/Colvars/spiff /tmp/spiff && \
    make -C /tmp/spiff && \
    install /tmp/spiff/spiff /usr/local/bin/ && \
    rm -fr /tmp/spiff

# TODO: there is no pre-built libtorch for arm64

# Install GitHub CLI
RUN dnf -y config-manager --add-repo https://cli.github.com/packages/rpm/gh-cli.repo && \
    dnf -y install gh --repo gh-cli && \
    dnf clean all

# Download pre-built recent CMake (mostly for GROMACS)
RUN mkdir -p /opt/cmake && \
    wget -O /tmp/cmake-install.sh https://github.com/Kitware/CMake/releases/download/v3.31.5/cmake-3.31.5-linux-aarch64.sh && \
    bash /tmp/cmake-install.sh --skip-license --prefix=/opt/cmake && \
    rm -f /tmp/cmake-install.sh

CMD ["/bin/bash", "-l"]
