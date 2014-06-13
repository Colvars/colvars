# Makefile for colvars plugin

.SUFFIXES: .C .o .d
.PHONY: clean depends

AR= ar
ARFLAGS = cr
RANLIB = ranlib

PLUGINNAME=colvars
VERSION=latest.version

COMPILEDIR = ../compile
ARCHDIR=$(COMPILEDIR)/lib_$(ARCH)/tcl/$(PLUGINNAME)
BINFILES=$(ARCHDIR)/colvars.so
TCLFILES=pkgIndex.tcl colvars.tcl 

SRCDIR=src
INCDIR=-Isrc -I$(VMDSRC) -I../include

VPATH = src $(ARCHDIR)

##
## Only build if we have a Tcl library 
##
ifdef TCLLIB
ifdef TCLINC
ifdef TCLLDFLAGS
TARGETS = $(ARCHDIR) $(BINFILES)
endif
endif
endif


# # # uncomment below for testing purposesARCH = LINUXAMD64
# VMDSRC = $(HOME)/programs/vmd/cvs/src
# ARCHDIR = obj
# TARGETS = $(ARCHDIR) $(BINFILES)
# CXX = g++
# CXXFLAGS = -m64 -O2 -fPIC
# COPTO = -fPIC -m64 -o 
# LOPTO = -fPIC -m64 -o 
# SHLD = gcc -shared
# # "ARCH = LINUXAMD64" \
# # "COPTO = -fPIC -m64 -o " \
# # "LOPTO = -fPIC -m64 -o " \
# # "CC = gcc" \
# # "CXX = g++" \
# # "DEF = -D" \
# # "TCLLDFLAGS = -ltcl8.5 -ldl" \
# # "NETCDFLDFLAGS = -lnetcdf " \
# # "AR = ar" \
# # "NM = nm -p" \
# # "RANLIB = touch" \
# # "SHLD = gcc -shared"


bins:
win32bins:
dynlibs: $(TARGETS) $(ARCHDIR)
staticlibs: 
win32staticlibs: 

distrib:
	for localname in `find ../compile -name colvars.so -print` ; do \
	pluginname=`echo $$localname | sed s/..\\\/compile\\\/lib_// `; \
	dir=`dirname $(PLUGINDIR)/$$pluginname`; \
	mkdir -p $$dir; \
	cp $$localname $(PLUGINDIR)/$$pluginname; \
	cp $(TCLFILES) $$dir ; \
	done

$(ARCHDIR):
	mkdir -p $(ARCHDIR)

COLVARS_SRCS = \
	src/colvaratoms.C \
	src/colvarbias.C \
	src/colvarbias_abf.C \
	src/colvarbias_meta.C \
	src/colvar.C \
	src/colvarcomp.C \
	src/colvarcomp_angles.C \
	src/colvarcomp_coordnums.C \
	src/colvarcomp_distances.C \
	src/colvarcomp_protein.C \
	src/colvarcomp_rotations.C \
	src/colvargrid.C \
	src/colvarmodule.C \
	src/colvarparse.C \
	src/colvarvalue.C \
	src/colvarproxy_vmd.C

## not ready yet
# COLVARS_SRCS + = src/colvarscript.C \
# 	src/colvarscript_tcl.C \

COLVARS_OBJS = $(patsubst %.C, %.o,  $(patsubst src/%,$(ARCHDIR)/%,$(COLVARS_SRCS)))

COLVARS_DEPS = $(COLVARS_SRCS:.C=.d)
.C.d:
	$(CXX) $(CXXFLAGS) $(INCDIR) -MM $< -MT $(patsubst src/%.C,'$$(ARCHDIR)'/%.o,$<) -MF $@

depends: Makefile.colvars 

Makefile.colvars: $(COLVARS_DEPS) Makefile 
	rm -f Makefile.colvars
	echo "COLVARS_OBJS = $(COLVARS_OBJS)" > Makefile.colvars
	for dep_file in src/*.d ; do \
	echo >> Makefile.colvars ; \
	cat $$dep_file >> Makefile.colvars ; \
	echo '	$$(CXX) $$(CXXFLAGS) $$(TCLINC) $$(INCDIR) -c $$< $$(COPTO)$$@' >> Makefile.colvars ;\
	done

# includes an explicit definition of COLVARS_OBJS
include Makefile.colvars

clean:
	rm -f ${COLVARS_OBJS} ${BINFILES} $(ARCHDIR)/colvars.so

$(ARCHDIR)/colvars.so: $(COLVARS_OBJS)
	if [ -n "$(TCLSHLD)" ]; \
	then $(TCLSHLD) $(LOPTO)$@ $(COLVARS_OBJS) $(TCLLIB) $(TCLLDFLAGS) $(LDFLAGS) $(CXXLDFLAGS) ; \
	else $(SHLD) $(LOPTO)$@ $(COLVARS_OBJS) $(TCLLIB) $(TCLLDFLAGS) $(CXXLDFLAGS) $(LDFLAGS); \
	fi


