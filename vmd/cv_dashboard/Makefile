
VERSION=1.0
PLUGINDIR=${HOME}/lib/vmd/plugins/noarch/tcl

install:
	mkdir -p ${PLUGINDIR}/cv_dashboard${VERSION} && \
		rsync -avL --delete cv_dashboard.tcl pkgIndex.tcl templates ${PLUGINDIR}/cv_dashboard${VERSION}/. 
