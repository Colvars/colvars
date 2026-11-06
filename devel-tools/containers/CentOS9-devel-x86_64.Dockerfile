FROM quay.io/centos/centos:stream9

LABEL org.opencontainers.image.title="CentOS Stream 9 development environment (x86_64 architecture)"

# Ensure Bash is used in the commands below
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
        gcc gcc-c++ gcc-gfortran \
        gcc-toolset-{12,13,14,15} \
        clang{,19,20}{,-devel,-analyzer,-tools-extra} \
        cppcheck libasan \
        autoconf automake cvs git subversion mercurial \
        cmake ccache ninja-build mold \
        Lmod openmpi-devel tbb-devel fftw-devel tcl-devel \
        python3-{devel,tkinter,virtualenv,numpy,scipy} \
        ncurses-devel \
        rlwrap \
        tk-devel fltk-devel \
        mesa-dri-drivers mesa-libGL-devel libglvnd-devel \
        libXi-devel libXcursor-devel libXinerama-devel libpng-devel \
        netcdf-devel hdf5-devel sqlite-devel expat-devel \
        doxygen \
        texlive-latex texlive-epstopdf texlive-pdftex texlive-tex4ht \
        texlive-graphics texlive-metafont texlive-collection-latexrecommended \
        texlive-collection-fontsrecommended \
        latexmk latex2html \
        inkscape gimp wget unzip \
        squashfs-tools && \
    dnf clean all

# The Clang package installs environment-modules, which conflicts with Lmod
RUN alternatives --set modules.sh /usr/share/lmod/lmod/init/profile

# Install GitHub CLI
RUN dnf -y config-manager --add-repo https://cli.github.com/packages/rpm/gh-cli.repo && dnf -y install gh --repo gh-cli

# Build Charm++
RUN umask 022 && \
    source /etc/profile && \
    module load mpi && \
    git clone --single-branch --depth=1 -b v8.0.1 https://github.com/UIUC-PPL/charm.git /opt/charm && \
    cd /opt/charm && \
    export CCACHE_DIR=/tmp && \
    ./build charm++ mpi-linux-x86_64 smp -j16 --with-production && \
    ./build charm++ multicore-linux-x86_64 -j16 --with-production && \
    ./build charm++ netlrts-linux-x86_64 smp -j16 --with-production

# Build and install spiff
RUN git clone https://github.com/Colvars/spiff /tmp/spiff && \
    make -C /tmp/spiff && \
    install /tmp/spiff/spiff /usr/local/bin/ && \
    rm -fr /tmp/spiff

# Download pre-built libTorch
RUN curl -L -o /tmp/libtorch.zip https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-2.4.1%2Bcpu.zip && \
    unzip /tmp/libtorch.zip -d /opt && \
    rm -f /tmp/libtorch.zip

# Download pre-built executable of recent CMake
RUN mkdir -p /opt/cmake && \
    curl -L -o /tmp/cmake-install.sh https://github.com/Kitware/CMake/releases/download/v4.4.5/cmake-4.4.5-linux-x86_64.sh && \
    bash /tmp/cmake-install.sh --skip-license --prefix=/opt/cmake && \
    rm -f /tmp/cmake-install.sh

CMD ["/bin/bash", "-l"]
