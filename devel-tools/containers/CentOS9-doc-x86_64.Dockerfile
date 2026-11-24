FROM quay.io/centos/centos:stream9

LABEL org.opencontainers.image.title="CentOS Stream 9 environment to build the Colvars documentation (x86_64 architecture)"

# Ensure Bash is used in the commands below
SHELL ["/bin/bash", "-c"]

# Same repositories as CentOS9-devel-x86_64, but only the tools needed by
# doc/Makefile: PDF (latexmk), HTML (LaTeXML), images (ImageMagick and
# Ghostscript), Doxygen (with Graphviz).  LaTeXML's compiled Perl dependencies
# come from the distribution; LaTeXML itself and TeX Live are installed below.
RUN dnf -y update && \
    dnf -y install dnf-plugins-core && \
    dnf config-manager --set-enabled crb && \
    dnf -y install \
        https://dl.fedoraproject.org/pub/epel/epel-release-latest-9.noarch.rpm \
        https://dl.fedoraproject.org/pub/epel/epel-next-release-latest-9.noarch.rpm && \
    dnf -y install \
        git make which diffutils findutils rsync wget unzip xz fontconfig \
        python3 \
        ImageMagick ghostscript \
        doxygen graphviz \
        perl perl-App-cpanminus perl-ExtUtils-MakeMaker perl-Test-Simple \
        perl-Archive-Zip perl-DB_File perl-File-Which perl-Getopt-Long \
        perl-JSON-XS perl-JSON-PP perl-libwww-perl perl-LWP-Protocol-https \
        perl-MIME-Base64 \
        perl-Parse-RecDescent perl-Pod-Parser perl-Text-Unidecode \
        perl-Time-HiRes perl-URI perl-XML-LibXML perl-XML-LibXSLT \
        ImageMagick-perl && \
    dnf clean all

# TeX Live from upstream: the distribution's TeX Live 2020 lacks packages used
# by the manual.  The frozen final snapshot of a given year makes the image
# reproducible; to upgrade, change TEXLIVE_YEAR.
ARG TEXLIVE_YEAR=2025
ARG TEXLIVE_REPO=https://ftp.math.utah.edu/pub/tex/historic/systems/texlive/${TEXLIVE_YEAR}/tlnet-final
ENV PATH=/opt/texlive/${TEXLIVE_YEAR}/bin/x86_64-linux:${PATH}

# Curated list of packages: those loaded by the manual (see the \usepackage
# lines in doc/colvars-refman*.tex) and their dependencies, plus latexmk.
# Check with "make clean-all && make" in doc/ after changing it.
RUN mkdir -p /tmp/install-tl && \
    curl -fsSL ${TEXLIVE_REPO}/install-tl-unx.tar.gz | \
        tar -xz -C /tmp/install-tl --strip-components=1 && \
    printf '%s\n' \
        "selected_scheme scheme-basic" \
        "TEXDIR /opt/texlive/${TEXLIVE_YEAR}" \
        "TEXMFLOCAL /opt/texlive/texmf-local" \
        "TEXMFSYSCONFIG /opt/texlive/${TEXLIVE_YEAR}/texmf-config" \
        "TEXMFSYSVAR /opt/texlive/${TEXLIVE_YEAR}/texmf-var" \
        "binary_x86_64-linux 1" \
        "instopt_adjustpath 0" \
        "tlpdbopt_autobackup 0" \
        "tlpdbopt_install_docfiles 0" \
        "tlpdbopt_install_srcfiles 0" \
        > /tmp/install-tl/texlive.profile && \
    /tmp/install-tl/install-tl --repository ${TEXLIVE_REPO} \
        --profile /tmp/install-tl/texlive.profile && \
    rm -rf /tmp/install-tl && \
    tlmgr option repository ${TEXLIVE_REPO} && \
    tlmgr install \
        latexmk \
        collection-latexrecommended collection-fontsrecommended \
        babel-english bbm bbm-macros catchfile comment etoolbox footmisc \
        jknapltx lipsum mdframed multirow natbib needspace rsfs \
        transparent zref && \
    chmod -R a+rX /opt/texlive

# LaTeXML (pure Perl; its compiled dependencies are installed above)
ARG LATEXML_VERSION=0.8.8
RUN (cpanm --notest --no-man-pages LaTeXML@${LATEXML_VERSION} || \
        (tail -50 ~/.cpanm/build.log; false)) && \
    rm -rf ~/.cpanm

# The CI jobs run as a different user, whose home may not be writable
ENV TEXMFVAR=/tmp/texmf-var

# Sanity checks
RUN latexmlc --VERSION && \
    latexmk -v | head -1 && \
    kpsewhich mdframed.sty bbm.sty mathrsfs.sty natbib.sty transparent.sty \
        mathptmx.sty helvet.sty && \
    convert -version | head -1 && \
    doxygen --version && \
    dot -V

CMD ["/bin/bash", "-l"]
