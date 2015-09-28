# pass version/platform information to compile
NAMD_VERSION = 2.10

# compiler flags (Win32 overrides)
COPTI = -I
COPTC = -c
COPTD = -D
COPTO = -o $(SPACE)

# Unix commands

ECHO = echo
MOVE = mv
MKDIR = mkdir -p
COPY = cp
RM = rm -f
LDD = ldd

# pass version/platform information to compile
RELEASE=$(COPTD)NAMD_VERSION=\"$(NAMD_VERSION)\" $(COPTD)NAMD_PLATFORM=\"$(NAMD_PLATFORM)\" $(RELEASEFLAGS)

# directories
SRCDIR = src
DSTDIR = obj
INCDIR = inc
DPMTADIR=dpmta-2.6
DPMEDIR=dpme2
PLUGINSRCDIR= plugins/molfile_plugin/src
PLUGININCDIR= plugins/include
SBSRCDIR = sb/src

MKDSTDIR = $(DSTDIR)/.exists
MKINCDIR = $(INCDIR)/.exists

# comment/uncomment these lines for (D)PMTA routines
#DPMTAINCL=$(COPTI)$(DPMTADIR)/mpole $(COPTI)$(DPMTADIR)/src
#DPMTALIB=-L$(DPMTADIR)/mpole -L$(DPMTADIR)/src -ldpmta2 -lmpole -lpvmc
#DPMTAFLAGS=$(COPTD)DPMTA
#DPMTA=$(DPMTAINCL) $(DPMTAFLAGS)
#DPMTALIBS=$(DPMTADIR)/mpole/libmpole.a $(DPMTADIR)/src/libdpmta2.a

# comment/uncomment these lines for DPME routines
#DPMEINCL=$(COPTI)$(DPMEDIR)
#DPMELIB=-L$(DPMEDIR) -ldpme
#DPMEFLAGS=$(COPTD)DPME
#DPME=$(DPMEINCL) $(DPMEFLAGS)
#DPMELIBS= $(DPMEDIR)/libdpme.a

# comment/uncomment these lines for FMM routines
#
# ufmmlap library from J. Huang - http://fastmultipole.org/Main/FMMSuite/
#   (ufmmlap = Uniform FMM Laplace Solver)
#
# Options below assume building with Gnu compilers.
# Define FMMDIR in Make.config file.
#
#FMMNAME=ufmmlap
#FMMINCL=$(COPTI)$(FMMDIR)/src
#FMMLIB=-L$(FMMDIR)/src -l$(FMMNAME) -lgfortran
#FMMFLAGS=$(COPTD)FMM_SOLVER
#FMM=$(FMMINCL) $(FMMFLAGS)
#FMMLIBS=$(FMMDIR)/src/lib$(FMMNAME).a

# to compile a memory optimized version, uncomment or config --with-memopt
#MEMOPT=-DMEM_OPT_VERSION
# to compile version that uses node-aware proxy send/recv spanning tree,
# add -DNODEAWARE_PROXY_SPANNINGTREE to the variable EXTRADEFINES
#EXTRADEFINES=-DREMOVE_PROXYDATAMSG_EXTRACOPY -DREMOVE_PROXYRESULTMSG_EXTRACOPY
EXTRADEFINES=-DREMOVE_PROXYRESULTMSG_EXTRACOPY -DNODEAWARE_PROXY_SPANNINGTREE -DUSE_NODEPATCHMGR
EXTRAINCS=
EXTRALINKLIBS=
# to compile namd using PAPI counters to measure flops and modify include and library path
# correspondingly
#EXTRADEFINES=-DREMOVE_PROXYRESULTMSG_EXTRACOPY -DMEASURE_NAMD_WITH_PAPI
#EXTRAINCS=-I$(HOME)/papi/include
#EXTRALINKLIBS=-lpapi -L$(HOME)/papi/lib

#If using the CkLoop library from Charm++ for PME calculation, first define macro
# USE_CKLOOP=1 for compiling the code, and then add "-module CkLoop" for linking.
CKLOOP = -DUSE_CKLOOP=1
CKLOOP_MODULE = -module CkLoop

# defaults for special cases
CXXTHREADOPTS = $(CXXOPTS) 
CXXSIMPARAMOPTS = $(CXXOPTS) 
CXXNOALIASOPTS = $(CXXOPTS) 
CXXMEMUSAGE = $(CXX)
CUDACC = $(CXX)
CUDAOBJS =
NATIVEPATH = echo

include Make.config

# define below Make.config so Win32 can change default target to winall
default: all

# Add new source files here.

OBJS = \
	$(DSTDIR)/ComputeMoa.o \
	$(DSTDIR)/ComputeMsm.o \
	$(DSTDIR)/ComputeMsmMsa.o \
	$(DSTDIR)/ComputeMsmSerial.o \
        $(DSTDIR)/ComputeFmmSerial.o \
	$(DSTDIR)/msm.o \
	$(DSTDIR)/msm_longrng.o \
	$(DSTDIR)/msm_longrng_sprec.o \
	$(DSTDIR)/msm_setup.o \
	$(DSTDIR)/msm_shortrng.o \
	$(DSTDIR)/msm_shortrng_sprec.o \
	$(DSTDIR)/wkfutils.o \
	$(DSTDIR)/common.o \
	$(DSTDIR)/dcdlib.o \
	$(DSTDIR)/erf.o \
	$(DSTDIR)/fitrms.o \
	$(DSTDIR)/main.o \
	$(DSTDIR)/mainfunc.o \
	$(DSTDIR)/memusage.o \
	$(DSTDIR)/strlib.o \
	$(DSTDIR)/AlgSeven.o \
	$(DSTDIR)/AlgRecBisection.o \
	$(DSTDIR)/AlgNbor.o \
	$(DSTDIR)/AtomMap.o \
	$(DSTDIR)/BackEnd.o \
	$(DSTDIR)/BroadcastMgr.o \
	$(DSTDIR)/BroadcastClient.o \
	$(DSTDIR)/CollectionMaster.o \
	$(DSTDIR)/CollectionMgr.o \
	$(DSTDIR)/Communicate.o \
	$(DSTDIR)/Compute.o \
	$(DSTDIR)/ComputeAngles.o \
	$(DSTDIR)/ComputeAniso.o \
	$(DSTDIR)/ComputeBonds.o \
	$(DSTDIR)/ComputeConsForce.o \
	$(DSTDIR)/ComputeConsForceMsgs.o \
	$(DSTDIR)/ComputeCrossterms.o \
	$(DSTDIR)/ComputeCylindricalBC.o \
	$(DSTDIR)/ComputeDihedrals.o \
	$(DSTDIR)/ComputeDPME.o \
	$(DSTDIR)/ComputeDPMEMsgs.o \
	$(DSTDIR)/ComputeDPMTA.o \
	$(DSTDIR)/ComputeEField.o \
	$(DSTDIR)/ComputeEwald.o \
	$(DSTDIR)/ComputeExt.o \
	$(DSTDIR)/ComputeGBISser.o \
	$(DSTDIR)/ComputeGBIS.o \
	$(DSTDIR)/ComputeGromacsPair.o \
	$(DSTDIR)/ComputeLCPO.o \
	$(DSTDIR)/ComputeFullDirect.o \
	$(DSTDIR)/ComputeHomePatch.o \
	$(DSTDIR)/ComputeHomePatches.o \
	$(DSTDIR)/ComputeImpropers.o \
	$(DSTDIR)/ComputeGlobal.o \
	$(DSTDIR)/ComputeGlobalMsgs.o \
	$(DSTDIR)/ComputeGridForce.o \
	$(DSTDIR)/ComputeMap.o \
	$(DSTDIR)/ComputeMgr.o \
	$(DSTDIR)/ComputeNonbondedSelf.o \
	$(DSTDIR)/ComputeNonbondedPair.o \
	$(DSTDIR)/ComputeNonbondedUtil.o \
	$(DSTDIR)/ComputeNonbondedStd.o \
	$(DSTDIR)/ComputeNonbondedFEP.o \
	$(DSTDIR)/ComputeNonbondedGo.o \
	$(DSTDIR)/ComputeNonbondedTI.o \
	$(DSTDIR)/ComputeNonbondedLES.o \
	$(DSTDIR)/ComputeNonbondedPProf.o \
	$(DSTDIR)/ComputeNonbondedTabEnergies.o \
	$(DSTDIR)/ComputeNonbondedCUDA.o \
	$(DSTDIR)/ComputeNonbondedCUDAExcl.o \
	$(DSTDIR)/ComputeNonbondedMIC.o \
	$(DSTDIR)/ComputeNonbondedMICKernel.o \
	$(DSTDIR)/ComputePatch.o \
	$(DSTDIR)/ComputePatchPair.o \
	$(DSTDIR)/ComputePme.o \
	$(DSTDIR)/OptPme.o \
	$(DSTDIR)/OptPmeRealSpace.o \
	$(DSTDIR)/ComputeRestraints.o \
	$(DSTDIR)/ComputeSphericalBC.o \
	$(DSTDIR)/ComputeStir.o \
	$(DSTDIR)/ComputeTclBC.o \
	$(DSTDIR)/ComputeThole.o \
	$(DSTDIR)/ConfigList.o \
	$(DSTDIR)/Controller.o \
	$(DSTDIR)/CudaUtils.o \
	$(DSTDIR)/ccsinterface.o \
	$(DSTDIR)/DataStream.o \
	$(DSTDIR)/DeviceCUDA.o \
	$(DSTDIR)/DumpBench.o \
        $(DSTDIR)/FreeEnergyAssert.o \
        $(DSTDIR)/FreeEnergyGroup.o \
        $(DSTDIR)/FreeEnergyLambda.o \
        $(DSTDIR)/FreeEnergyLambdMgr.o \
        $(DSTDIR)/FreeEnergyParse.o \
        $(DSTDIR)/FreeEnergyRestrain.o \
        $(DSTDIR)/FreeEnergyRMgr.o \
        $(DSTDIR)/FreeEnergyVector.o \
	$(DSTDIR)/fstream_namd.o \
	$(DSTDIR)/GlobalMaster.o \
	$(DSTDIR)/GlobalMasterServer.o \
	$(DSTDIR)/GlobalMasterTest.o \
	$(DSTDIR)/GlobalMasterIMD.o \
	$(DSTDIR)/GlobalMasterTcl.o \
	$(DSTDIR)/GlobalMasterSMD.o \
	$(DSTDIR)/GlobalMasterTMD.o \
  $(DSTDIR)/Matrix4Symmetry.o \
  $(DSTDIR)/GlobalMasterSymmetry.o \
	$(DSTDIR)/GlobalMasterFreeEnergy.o \
	$(DSTDIR)/GlobalMasterEasy.o \
	$(DSTDIR)/GlobalMasterMisc.o \
	$(DSTDIR)/colvar.o \
	$(DSTDIR)/colvaratoms.o \
	$(DSTDIR)/colvarbias.o \
	$(DSTDIR)/colvarbias_abf.o \
	$(DSTDIR)/colvarbias_alb.o \
	$(DSTDIR)/colvarbias_histogram.o \
	$(DSTDIR)/colvarbias_meta.o \
	$(DSTDIR)/colvarbias_restraint.o \
	$(DSTDIR)/colvarcomp.o \
	$(DSTDIR)/colvarcomp_angles.o \
	$(DSTDIR)/colvarcomp_coordnums.o \
	$(DSTDIR)/colvarcomp_distances.o \
	$(DSTDIR)/colvarcomp_protein.o \
	$(DSTDIR)/colvarcomp_rotations.o \
	$(DSTDIR)/colvargrid.o \
	$(DSTDIR)/colvarmodule.o \
	$(DSTDIR)/colvarparse.o \
	$(DSTDIR)/colvarproxy_namd.o \
	$(DSTDIR)/colvarscript.o \
	$(DSTDIR)/colvartypes.o \
	$(DSTDIR)/colvarvalue.o \
	$(DSTDIR)/GridForceGrid.o \
        $(DSTDIR)/GromacsTopFile.o \
	$(DSTDIR)/heap.o \
	$(DSTDIR)/HomePatch.o \
	$(DSTDIR)/IMDOutput.o \
	$(DSTDIR)/InfoStream.o \
	$(DSTDIR)/LdbCoordinator.o \
	$(DSTDIR)/LJTable.o \
	$(DSTDIR)/Measure.o \
	$(DSTDIR)/MGridforceParams.o \
	$(DSTDIR)/MStream.o \
	$(DSTDIR)/MigrateAtomsMsg.o \
	$(DSTDIR)/Molecule.o \
	$(DSTDIR)/Molecule2.o \
	$(DSTDIR)/GoMolecule.o \
	$(DSTDIR)/NamdCentLB.o \
	$(DSTDIR)/NamdNborLB.o \
	$(DSTDIR)/NamdHybridLB.o \
	$(DSTDIR)/NamdDummyLB.o \
	$(DSTDIR)/NamdState.o \
	$(DSTDIR)/NamdOneTools.o \
	$(DSTDIR)/Node.o \
	$(DSTDIR)/Output.o \
	$(DSTDIR)/Parameters.o \
	$(DSTDIR)/ParseOptions.o \
	$(DSTDIR)/Patch.o \
	$(DSTDIR)/PatchMgr.o \
	$(DSTDIR)/PatchMap.o \
	$(DSTDIR)/PDB.o \
	$(DSTDIR)/PDBData.o \
	$(DSTDIR)/PmeKSpace.o \
	$(DSTDIR)/PmeRealSpace.o \
	$(DSTDIR)/ProcessorPrivate.o \
	$(DSTDIR)/ProxyMgr.o \
	$(DSTDIR)/ProxyPatch.o \
	$(DSTDIR)/Rebalancer.o \
	$(DSTDIR)/RecBisection.o \
	$(DSTDIR)/ReductionMgr.o \
	$(DSTDIR)/RefineOnly.o \
	$(DSTDIR)/RefineTorusLB.o \
	$(DSTDIR)/ScriptTcl.o \
	$(DSTDIR)/Sequencer.o \
	$(DSTDIR)/Set.o \
	$(DSTDIR)/Settle.o \
	$(DSTDIR)/SimParameters.o \
	$(DSTDIR)/SortAtoms.o \
	$(DSTDIR)/Sync.o \
	$(DSTDIR)/TclCommands.o \
	$(DSTDIR)/TorusLB.o \
	$(DSTDIR)/WorkDistrib.o \
	$(DSTDIR)/pub3dfft.o \
	$(DSTDIR)/vmdsock.o \
	$(DSTDIR)/parm.o \
	$(DSTDIR)/imd.o \
	$(DSTDIR)/CompressPsf.o \
	$(DSTDIR)/PluginIOMgr.o \
	$(DSTDIR)/DataExchanger.o \
	$(DSTDIR)/ParallelIOMgr.o 


# Add new modules here.

CIFILES = 	\
		$(INCDIR)/ComputeMoaMgr.decl.h \
		$(INCDIR)/ComputeMoaMgr.def.h \
		$(INCDIR)/ComputeMsmMgr.decl.h \
		$(INCDIR)/ComputeMsmMgr.def.h \
		$(INCDIR)/ComputeMsmMsaMgr.decl.h \
		$(INCDIR)/ComputeMsmMsaMgr.def.h \
		$(INCDIR)/ComputeMsmSerialMgr.decl.h \
		$(INCDIR)/ComputeMsmSerialMgr.def.h \
		$(INCDIR)/ComputeFmmSerialMgr.decl.h \
		$(INCDIR)/ComputeFmmSerialMgr.def.h \
		$(INCDIR)/BroadcastMgr.decl.h \
		$(INCDIR)/BroadcastMgr.def.h \
		$(INCDIR)/CollectionMaster.decl.h \
		$(INCDIR)/CollectionMaster.def.h \
		$(INCDIR)/CollectionMgr.decl.h \
		$(INCDIR)/CollectionMgr.def.h \
		$(INCDIR)/ComputeMgr.decl.h \
		$(INCDIR)/ComputeMgr.def.h \
		$(INCDIR)/ComputeGridForceMgr.decl.h \
		$(INCDIR)/ComputeGridForceMgr.def.h \
		$(INCDIR)/ComputePmeMgr.decl.h \
		$(INCDIR)/ComputePmeMgr.def.h \
		$(INCDIR)/OptPmeMgr.decl.h \
		$(INCDIR)/OptPmeMgr.def.h \
		$(INCDIR)/PmeFFTLib.decl.h \
		$(INCDIR)/PmeFFTLib.def.h \
		$(INCDIR)/ComputeExtMgr.decl.h \
		$(INCDIR)/ComputeExtMgr.def.h \
		$(INCDIR)/ComputeGBISserMgr.decl.h \
		$(INCDIR)/ComputeGBISserMgr.def.h \
		$(INCDIR)/LdbCoordinator.decl.h \
		$(INCDIR)/LdbCoordinator.def.h \
		$(INCDIR)/NamdCentLB.decl.h \
		$(INCDIR)/NamdCentLB.def.h \
		$(INCDIR)/NamdNborLB.decl.h \
		$(INCDIR)/NamdNborLB.def.h \
		$(INCDIR)/NamdHybridLB.decl.h \
		$(INCDIR)/NamdHybridLB.def.h \
		$(INCDIR)/NamdDummyLB.decl.h \
		$(INCDIR)/NamdDummyLB.def.h \
		$(INCDIR)/Node.decl.h \
		$(INCDIR)/Node.def.h \
		$(INCDIR)/PatchMgr.decl.h \
		$(INCDIR)/PatchMgr.def.h \
		$(INCDIR)/ProxyMgr.decl.h \
		$(INCDIR)/ProxyMgr.def.h \
		$(INCDIR)/ReductionMgr.decl.h \
		$(INCDIR)/ReductionMgr.def.h \
		$(INCDIR)/Sync.decl.h \
		$(INCDIR)/Sync.def.h \
		$(INCDIR)/WorkDistrib.decl.h \
		$(INCDIR)/WorkDistrib.def.h \
		$(INCDIR)/ParallelIOMgr.decl.h \
		$(INCDIR)/ParallelIOMgr.def.h \
		$(INCDIR)/DataExchanger.decl.h \
		$(INCDIR)/DataExchanger.def.h \
		$(INCDIR)/main.decl.h \
		$(INCDIR)/main.def.h 

# Add new source files here.

PLUGINOBJS = \
	$(DSTDIR)/dcdplugin.o \
	$(DSTDIR)/jsplugin.o \
	$(DSTDIR)/namdbinplugin.o \
	$(DSTDIR)/pdbplugin.o \
	$(DSTDIR)/psfplugin.o

PLUGINLIB = $(PLUGINOBJS)

CUDAOBJSRAW = \
	$(DSTDIR)/ComputePmeCUDAKernel.o \
	$(DSTDIR)/ComputeNonbondedCUDAKernel.o

$(DSTDIR)/ComputePmeCUDAKernel.o: \
	$(SRCDIR)/ComputePmeCUDAKernel.cu \
	$(SRCDIR)/ComputePmeCUDAKernel.h
	$(CUDACC) $(CUDACCOPTS) -Xptxas -v $(COPTO) "`$(NATIVEPATH) $(DSTDIR)/`ComputePmeCUDAKernel.o" $(COPTC) "`$(NATIVEPATH) $(SRCDIR)/`ComputePmeCUDAKernel.cu"

$(DSTDIR)/ComputeNonbondedCUDAKernel.o: \
	$(SRCDIR)/ComputeNonbondedCUDAKernel.cu \
	$(SRCDIR)/ComputeNonbondedCUDAKernelBase.h \
	$(SRCDIR)/ComputeGBISCUDAKernel.h \
	$(SRCDIR)/ComputeNonbondedCUDAKernel.h
	$(CUDACC) $(CUDACCOPTS) -Xptxas -v $(COPTO) "`$(NATIVEPATH) $(DSTDIR)/`ComputeNonbondedCUDAKernel.o" $(COPTC) "`$(NATIVEPATH) $(SRCDIR)/`ComputeNonbondedCUDAKernel.cu"

SBOBJS = \
	$(DSTDIR)/tcl_main.o \
	$(DSTDIR)/tcl_psfgen.o \
	$(DSTDIR)/charmm_file.o \
	$(DSTDIR)/charmm_parse_topo_defs.o \
	$(DSTDIR)/extract_alias.o \
	$(DSTDIR)/hash.o \
	$(DSTDIR)/hasharray.o \
	$(DSTDIR)/memarena.o \
	$(DSTDIR)/pdb_file.o \
	$(DSTDIR)/pdb_file_extract.o \
	$(DSTDIR)/psf_file.o \
	$(DSTDIR)/psf_file_extract.o \
	$(DSTDIR)/topo_defs.o \
	$(DSTDIR)/topo_mol.o \
	$(DSTDIR)/topo_mol_output.o \
	$(DSTDIR)/topo_mol_pluginio.o \
	$(DSTDIR)/stringhash.o

# definitions for Charm routines
CHARMC = $(CHARM)/bin/charmc
CHARMXI = $(CHARM)/bin/charmc
CHARMINC = $(CHARM)/include $(COPTD)CMK_OPTIMIZE=1
CHARMLIB = $(CHARM)/lib
CHARM_MODULES = -module NeighborLB -module HybridLB -module RefineLB -module GreedyLB -module CkMulticast $(CKLOOP_MODULE)
#CHARM_MODULES = -module NeighborLB -module HybridLB -module RefineLB -module GreedyLB
#CHARM_MODULES = -module msa -module NeighborLB -module HybridLB -module RefineLB -module GreedyLB
#MSA = -DCHARM_HAS_MSA

# Libraries we may have changed
LIBS = $(CUDAOBJS) $(PLUGINLIB) $(DPMTALIBS) $(DPMELIBS) $(FMMLIBS) $(TCLDLL)

# CXX is platform dependent
CXXBASEFLAGS = $(COPTI)$(CHARMINC) $(COPTI)$(SRCDIR) $(COPTI)$(INCDIR) $(DPMTA) $(DPME) $(FMM) $(COPTI)$(PLUGININCDIR) $(COPTD)STATIC_PLUGIN $(TCL) $(FFT) $(CUDA) $(MIC) $(MEMOPT) $(CCS) $(RELEASE) $(EXTRADEFINES) $(TRACEOBJDEF) $(EXTRAINCS) $(MSA) $(CKLOOP)
CXXFLAGS = $(CXXBASEFLAGS) $(CXXOPTS) $(CXXMICOPTS)
CXXTHREADFLAGS = $(CXXBASEFLAGS) $(CXXTHREADOPTS) $(CXXMICOPTS)
CXXSIMPARAMFLAGS = $(CXXBASEFLAGS) $(CXXSIMPARAMOPTS) $(CXXMICOPTS)
CXXNOALIASFLAGS = $(CXXBASEFLAGS) $(CXXNOALIASOPTS) $(CXXMICOPTS)
GXXFLAGS = $(CXXBASEFLAGS) -DNO_STRSTREAM_H
CFLAGS = $(COPTI)$(SRCDIR) $(TCL) $(COPTS) $(RELEASE) $(EXTRADEFINES) $(TRACEOBJDEF)
PLUGINGCCFLAGS = $(COPTI)$(PLUGINSRCDIR) $(COPTI)$(PLUGININCDIR) $(COPTD)STATIC_PLUGIN
PLUGINCFLAGS = $(PLUGINGCCFLAGS) $(COPTS)
SBCFLAGS = $(COPTI)$(SBSRCDIR) $(COPTI)$(PLUGININCDIR) $(COPTD)STATIC_PLUGIN -DPSFGEN_USEPLUGINS $(TCL) $(COPTS) $(RELEASE) $(EXTRADEFINES) $(TRACEOBJDEF)
SBGCCFLAGS = $(COPTI)$(SBSRCDIR) $(COPTI)$(PLUGININCDIR) $(COPTD)STATIC_PLUGIN -DPSFGEN_USEPLUGINS $(TCL) $(RELEASE) $(EXTRADEFINES) $(TRACEOBJDEF)

# Add new executables here.

BINARIES = namd2 psfgen sortreplicas flipdcd flipbinpdb charmrun

# This should be rebuilt at every compile, but not on Win32.
BUILDINFO = $(DSTDIR)/buildinfo
MAKEBUILDINFO = \
	$(RM) $(BUILDINFO).C; \
	echo 'const char *namd_build_date = ' \"`date`\"\; > $(BUILDINFO).C; \
	echo 'const char *namd_build_user = ' \"$(USER)\"\; >> $(BUILDINFO).C; \
	echo 'const char *namd_build_machine = ' \"`hostname`\"\; >> $(BUILDINFO).C; \
	cat $(BUILDINFO).C; \
	$(CXX) $(CXXFLAGS) $(COPTO)$(BUILDINFO).o $(COPTC) $(BUILDINFO).C

all:	$(BINARIES) $(LIBCUDARTSO)

namd2:	$(MKINCDIR) $(MKDSTDIR) $(OBJS) $(LIBS)
	$(MAKEBUILDINFO)
	$(CHARMC) -verbose -ld++-option \
	'$(COPTI)$(CHARMINC) $(COPTI)$(INCDIR) $(COPTI)$(SRCDIR) $(CXXOPTS) $(CXXMICOPTS)' \
	$(CHARM_MODULES) -language charm++ \
	$(BUILDINFO).o \
	$(OBJS) \
	$(CUDAOBJS) \
	$(CUDALIB) \
	$(DPMTALIB) \
	$(DPMELIB) \
	$(FMMLIB) \
	$(TCLLIB) \
	$(FFTLIB) \
	$(PLUGINLIB) \
	$(CHARMOPTS) \
	$(EXTRALINKLIBS) \
	-lm -o namd2

charmrun: $(CHARM)/bin/charmrun # XXX
	$(COPY) $(CHARM)/bin/charmrun $@

$(LIBCUDARTSO):
	$(COPY) $(CUDASODIR)/$(LIBCUDARTSO) $@;

WINDOWSBINARIES = namd2.exe psfgen.exe sortreplicas.exe

winall: $(WINDOWSBINARIES) $(LIBCUDARTSO)

namd2.exe:  $(MKINCDIR) $(MKDSTDIR) $(OBJS) $(LIBS) $(TCLDLL)
	$(MAKEBUILDINFO)
	$(CHARMC) -verbose \
	$(CHARM_MODULES) -language charm++ \
	$(BUILDINFO).o \
	$(OBJS) \
	$(CUDAOBJS) \
	$(CUDALIB) \
	$(DPMTALIB) \
	$(DPMELIB) \
	$(FMMLIB) \
	$(TCLLIB) \
	$(FFTLIB) \
	$(PLUGINLIB) \
	$(CHARMOPTS) \
	-o namd2

charmd.exe:
	$(COPY) $(CHARM)/bin/charmd.exe charmd.exe

charmd_faceless.exe:
	$(COPY) $(CHARM)/bin/charmd_faceless.exe charmd_faceless.exe

charmrun.exe:
	$(COPY) $(CHARM)/bin/charmrun.exe charmrun.exe

psfgen:	$(MKDSTDIR) $(SBOBJS) $(PLUGINOBJS)
	$(CC) $(SBCFLAGS) -o psfgen $(SBOBJS) $(PLUGINOBJS) $(TCLLIB) $(TCLAPPLIB) -lm

psfgen.exe:	$(MKDSTDIR) $(SBOBJS) $(PLUGINOBJS) $(TCLDLL)
	$(CC) $(SBCFLAGS) -o psfgen $(SBOBJS) $(PLUGINOBJS) $(TCLLIB) $(TCLAPPLIB) -lm

sortreplicas:	$(MKDSTDIR) $(DSTDIR)/sortreplicas.o $(PLUGINOBJS)
	$(CC) $(SBCFLAGS) -o sortreplicas $(DSTDIR)/sortreplicas.o $(PLUGINOBJS) -lm

sortreplicas.exe:	$(MKDSTDIR) $(DSTDIR)/sortreplicas.o $(PLUGINOBJS)
	$(CC) $(SBCFLAGS) -o sortreplicas $(DSTDIR)/sortreplicas.o $(PLUGINOBJS) -lm

$(DSTDIR)/sortreplicas.o:	$(MKDSTDIR) $(SRCDIR)/sortreplicas.c
	$(CC) $(SBCFLAGS) $(COPTO)$(DSTDIR)/sortreplicas.o $(COPTC) $(SRCDIR)/sortreplicas.c

diffbinpdb:	$(SRCDIR)/diffbinpdb.c
	$(CC) $(CFLAGS) -o diffbinpdb $(SRCDIR)/diffbinpdb.c -lm

flipdcd:	$(SRCDIR)/flipdcd.c
	$(CC) $(CFLAGS) -o $@ $(SRCDIR)/flipdcd.c || \
	echo "#!/bin/sh\necho unavailable on this platform" > $@; \
	chmod +x $@

flipbinpdb:	$(SRCDIR)/flipbinpdb.c
	$(CC) $(CFLAGS) -o $@ $(SRCDIR)/flipbinpdb.c || \
	echo "#!/bin/sh\necho unavailable on this platform" > $@; \
	chmod +x $@

fixdcd:	$(SRCDIR)/fixdcd.c
	$(CC) $(CFLAGS) -o fixdcd $(SRCDIR)/fixdcd.c

dumpdcd:	$(SRCDIR)/dumpdcd.c
	$(CC) $(CFLAGS) -o dumpdcd $(SRCDIR)/dumpdcd.c

loaddcd:	$(SRCDIR)/loaddcd.c
	$(CC) $(CFLAGS) -o loaddcd $(SRCDIR)/loaddcd.c

updatefiles:
	touch ../src/ComputeSelfTuples.h
	rm -f obj/ComputeNonbondedPair.o
	rm -f obj/ComputeNonbondedSelf.o
	rm -f obj/ComputePme.o
	rm -f obj/OptPme.o

#To compile tracecomputes, type the command "make tracecomputes TRACEOBJDEF=-DTRACE_COMPUTE_OBJECTS"
tracecomputes: updatefiles $(MKINCDIR) $(MKDSTDIR) $(OBJS) $(LIBS)
	$(MAKEBUILDINFO)
	$(CHARMC) -verbose -ld++-option \
	'$(COPTI)$(CHARMINC) $(COPTI)$(INCDIR) $(COPTI)$(SRCDIR) $(CXXOPTS) $(CXXMICOPTS)' \
	$(CHARM_MODULES) -language charm++ \
	-tracemode projections \
	$(BUILDINFO).o \
	$(OBJS) \
	$(CUDAOBJS) \
	$(CUDALIB) \
	$(DPMTALIB) \
	$(DPMELIB) \
	$(FMMLIB) \
	$(TCLLIB) \
	$(FFTLIB) \
	$(PLUGINLIB) \
	$(CHARMOPTS) \
	$(EXTRALINKLIBS) \
	-lm -o namd2.tc.prj

projections: $(MKINCDIR) $(MKDSTDIR) $(OBJS) $(LIBS)
	$(MAKEBUILDINFO)
	$(CHARMC) -verbose -ld++-option \
	'$(COPTI)$(CHARMINC) $(COPTI)$(INCDIR) $(COPTI)$(SRCDIR) $(CXXOPTS) $(CXXMICOPTS)' \
	$(CHARM_MODULES) -language charm++ \
	-tracemode projections \
	$(BUILDINFO).o \
	$(OBJS) \
	$(CUDAOBJS) \
	$(CUDALIB) \
	$(DPMTALIB) \
	$(DPMELIB) \
	$(FMMLIB) \
	$(TCLLIB) \
	$(FFTLIB) \
	$(PLUGINLIB) \
	$(CHARMOPTS) \
	$(EXTRALINKLIBS) \
	-lm -o namd2.prj

summary: $(MKINCDIR) $(MKDSTDIR) $(OBJS) $(LIBS)
	$(MAKEBUILDINFO)
	$(CHARMC) -verbose -ld++-option \
	'$(COPTI)$(CHARMINC) $(COPTI)$(INCDIR) $(COPTI)$(SRCDIR) $(CXXOPTS) $(CXXMICOPTS)' \
	-module NeighborLB -module CkMulticast -language charm++ \
	-tracemode summary \
	$(BUILDINFO).o \
	$(OBJS) \
	$(CUDAOBJS) \
	$(CUDALIB) \
	$(DPMTALIB) \
	$(DPMELIB) \
	$(FMMLIB) \
	$(TCLLIB) \
	$(FFTLIB) \
	$(PLUGINLIB) \
	$(CHARMOPTS) \
	$(EXTRALINKLIBS) \
	-lm -o namd2.sum

$(DPMTADIR)/mpole/libmpole.a: $(DPMTADIR)/src/libdpmta2.a

$(DPMTADIR)/src/libdpmta2.a:
	cd $(DPMTADIR) ; $(MAKE) ; cd ..

$(DPMEDIR)/libdpme.a:
	cd $(DPMEDIR) ; $(MAKE) ; cd ..

$(FMMDIR)/src/lib$(FMMNAME).a:
	cd $(FMMDIR) ; $(MAKE) ; cd ..


# Implicit rules for modules.

.SECONDARY:
	# prevent gmake from deleting intermediate files

$(INCDIR)/%.decl.h $(INCDIR)/%.def.h: $(MKINCDIR) $(SRCDIR)/%.ci
	$(COPY) $(SRCDIR)/$*.ci $(INCDIR)
	$(CHARMXI) $(INCDIR)/$*.ci
	$(RM) $(INCDIR)/$*.ci
	$(MOVE) $*.def.h $(INCDIR)
	$(MOVE) $*.decl.h $(INCDIR)


# Explicit rules for modules that don't match their file names.
# Multiple targets must be a pattern to execute recipe only once.

$(INCDIR)/PmeFF%Lib.decl.h $(INCDIR)/PmeFF%Lib.def.h: $(MKINCDIR) $(SRCDIR)/fftlib.ci
	$(COPY) $(SRCDIR)/fftlib.ci $(INCDIR)
	$(CHARMXI) $(INCDIR)/fftlib.ci
	$(RM) $(INCDIR)/fftlib.ci
	$(MOVE) PmeFFTLib.def.h $(INCDIR)
	$(MOVE) PmeFFTLib.decl.h $(INCDIR)

$(INCDIR)/OptPm%Mgr.decl.h $(INCDIR)/OptPm%Mgr.def.h: $(MKINCDIR) $(SRCDIR)/OptPme.ci
	$(COPY) $(SRCDIR)/OptPme.ci $(INCDIR)
	$(CHARMXI) $(INCDIR)/OptPme.ci
	$(RM) $(INCDIR)/OptPme.ci
	$(MOVE) OptPmeMgr.def.h $(INCDIR)
	$(MOVE) OptPmeMgr.decl.h $(INCDIR)

DEPENDFILE = .rootdir/Make.depends

# This is a CPU killer...  Don't make depends if you don't need to.
depends: $(MKINCDIR) $(CIFILES) $(MKDSTDIR) $(DEPENDFILE)
	$(ECHO) "Creating " $(DEPENDFILE) " ..."; \
	if [ -f $(DEPENDFILE) ]; then \
	   $(MOVE) -f $(DEPENDFILE) $(DEPENDFILE).old; \
	fi; \
	touch $(DEPENDFILE); \
	for i in $(OBJS) ; do \
	      SRCFILE=$(SRCDIR)/`basename $$i .o`.C ; \
	      if [ ! -f $$SRCFILE ]; then \
	            SRCFILE=$(SRCDIR)/`basename $$i .o`.c ; \
              fi; \
	      $(ECHO) "checking dependencies for $$SRCFILE" ; \
	      g++ -MM $(GXXFLAGS) $$SRCFILE | \
	      perl $(SRCDIR)/dc.pl $(CHARMINC) $(TCLDIR) $(FFTDIR) /usr/include /usr/local >> $(DEPENDFILE); \
	      $(ECHO) '	$$(CXX) $$(CXXFLAGS) $$(COPTO)'$$i '$$(COPTC)' \
		$$SRCFILE >> $(DEPENDFILE) ; \
	done; \
	for i in $(PLUGINOBJS) ; do \
	      BASENAME=`basename $$i .o` ; \
	      SRCFILE=$(PLUGINSRCDIR)/$$BASENAME.c ; \
	      $(ECHO) "checking dependencies for $$SRCFILE" ; \
	      gcc -MM $(PLUGINGCCFLAGS) $$SRCFILE | \
	      perl $(SRCDIR)/dc.pl /usr/include /usr/local >> $(DEPENDFILE); \
	      $(ECHO) '	$$(CC) $$(PLUGINCFLAGS) $$(COPTO)'$$i '$$(COPTC)' \
		'$$(COPTD)'VMDPLUGIN=molfile_$$BASENAME \
		$$SRCFILE >> $(DEPENDFILE) ; \
	done; \
	for i in $(SBOBJS) ; do \
	      SRCFILE=$(SBSRCDIR)/`basename $$i .o`.c ; \
	      $(ECHO) "checking dependencies for $$SRCFILE" ; \
	      gcc -MM $(SBGCCFLAGS) $$SRCFILE | \
	      perl $(SRCDIR)/dc.pl $(CHARMINC) $(TCLDIR) $(FFTDIR) /usr/include /usr/local >> $(DEPENDFILE); \
	      $(ECHO) '	$$(CC) $$(SBCFLAGS) $$(COPTO)'$$i '$$(COPTC)' \
		$$SRCFILE >> $(DEPENDFILE) ; \
	done; \
	$(RM) $(DEPENDFILE).sed; \
	sed -e "/obj\/Controller.o/ s/CXXFLAGS/CXXTHREADFLAGS/" \
	    -e "/obj\/Sequencer.o/ s/CXXFLAGS/CXXTHREADFLAGS/" \
	    -e "/obj\/ComputeFullDirect.o/ s/CXXFLAGS/CXXTHREADFLAGS/" \
	    -e "/obj\/ReductionMgr.o/ s/CXXFLAGS/CXXTHREADFLAGS/" \
	    -e "/obj\/SimParameters.o/ s/CXXFLAGS/CXXSIMPARAMFLAGS/" \
	    -e "/obj\/ComputeNonbondedStd.o/ s/CXXFLAGS/CXXNOALIASFLAGS/" \
	    -e "/obj\/ComputeNonbondedFEP.o/ s/CXXFLAGS/CXXNOALIASFLAGS/" \
	    -e "/obj\/ComputeNonbondedTI.o/ s/CXXFLAGS/CXXNOALIASFLAGS/" \
	    -e "/obj\/ComputeNonbondedLES.o/ s/CXXFLAGS/CXXNOALIASFLAGS/" \
	    -e "/obj\/ComputeNonbondedPProf.o/ s/CXXFLAGS/CXXNOALIASFLAGS/" \
	    -e "/obj\/ComputeNonbondedTabEnergies.o/ s/CXXFLAGS/CXXNOALIASFLAGS/" \
	    -e "/obj\/memusage.o/ s/CXX/CXXMEMUSAGE/" \
	    $(DEPENDFILE) > $(DEPENDFILE).sed; \
	$(MOVE) -f $(DEPENDFILE).sed $(DEPENDFILE);

$(DEPENDFILE):
	touch $(DEPENDFILE)

include Make.depends

$(MKDSTDIR):
	if [ ! -d $(DSTDIR) ]; then $(MKDIR) $(DSTDIR); fi
	if [ ! -f $(MKDSTDIR) ]; then touch $(MKDSTDIR); fi

$(MKINCDIR):
	if [ ! -d $(INCDIR) ]; then $(MKDIR) $(INCDIR); fi
	if [ ! -f $(MKINCDIR) ]; then touch $(MKINCDIR); fi

clean:
	rm -rf ptrepository Templates.DB SunWS_cache $(DSTDIR) $(INCDIR)

veryclean:	clean
	rm -f $(BINARIES)

RELEASE_DIR_NAME = NAMD_$(NAMD_VERSION)_$(NAMD_PLATFORM)

DOC_FILES = README.txt announce.txt license.txt notes.txt

RELEASE_FILES = $(LIBCUDARTSO) flipdcd flipbinpdb sortreplicas psfgen charmrun namd2

WINDOWS_RELEASE_FILES = $(LIBCUDARTSO) $(WINDOWSBINARIES) $(TCLDLL)

release: all
	$(ECHO) Creating release $(RELEASE_DIR_NAME)
	mkdir $(RELEASE_DIR_NAME)
	cp $(RELEASE_FILES) $(COMPILERRUNTIMELIBS) $(RELEASE_DIR_NAME)
	for f in $(DOC_FILES); do cp .rootdir/$$f $(RELEASE_DIR_NAME); done
	cp -r .rootdir/lib $(RELEASE_DIR_NAME)
	for f in `find $(RELEASE_DIR_NAME)/lib -name CVS`; do \
	  /bin/rm -rf $$f; \
	done
	if [ -r $(CHARM)/bin/charmd ]; then \
	  $(COPY) $(CHARM)/bin/charmd $(RELEASE_DIR_NAME); \
	fi
	if [ -r $(CHARM)/bin/charmd_faceless ]; then \
	  $(COPY) $(CHARM)/bin/charmd_faceless $(RELEASE_DIR_NAME); \
	fi
	chmod -R a+rX $(RELEASE_DIR_NAME)
	tar cf $(RELEASE_DIR_NAME).tar $(RELEASE_DIR_NAME)
	gzip $(RELEASE_DIR_NAME).tar
	echo $(CHARM)
	ls -l $(CHARM)/lib
	-for f in $(RELEASE_FILES); do echo $$f; $(LDD) $(RELEASE_DIR_NAME)/$$f; done

winrelease: winall
	$(ECHO) Creating release $(RELEASE_DIR_NAME)
	mkdir $(RELEASE_DIR_NAME)
	cp $(WINDOWS_RELEASE_FILES) $(RELEASE_DIR_NAME)
	for f in $(DOC_FILES); do \
	  perl -p -i -e 's/(?<!\r)\n$$/\r\n/' < .rootdir/$$f > $(RELEASE_DIR_NAME)/$$f; \
	done
	cp -r .rootdir/lib $(RELEASE_DIR_NAME)
	for f in `find $(RELEASE_DIR_NAME)/lib -name CVS`; do \
	  /bin/rm -rf $$f; \
	done
	for f in `find $(RELEASE_DIR_NAME)/lib -type f`; do \
	  perl -p -i -e 's/(?<!\r)\n$$/\r\n/' < $$f > $$f.wintxt; \
	  mv $$f.wintxt $$f; \
	done
	chmod -R a+rX $(RELEASE_DIR_NAME)
	echo $(CHARM)
	ls -l $(CHARM)/lib
	echo $(CHARM)
	zip -r $(RELEASE_DIR_NAME).zip $(RELEASE_DIR_NAME)

