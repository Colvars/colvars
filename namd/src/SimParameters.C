/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*****************************************************************************
 * $Source: /home/cvs/namd/cvsroot/namd2/src/SimParameters.C,v $
 * $Author: jim $
 * $Date: 2017/03/30 20:06:17 $
 * $Revision: 1.1478 $
 *****************************************************************************/

/** \file SimParameters.C
 * SimParameters is just a glorified structure to hold the global
 * static simulation parameters such as timestep size, cutoff, etc. that
 * are read in from the configuration file.
 **/

#include "InfoStream.h"
#include "ComputeNonbondedUtil.h"
#include "LJTable.h"
#include "ComputePme.h"
#include "ComputeCUDAMgr.h"
#include "ConfigList.h"
#include "SimParameters.h"
#include "ParseOptions.h"
#include "structures.h"
#include "Communicate.h"
#include "MStream.h"
#include "Output.h"
#include <stdio.h>
#include <time.h>
#ifdef NAMD_FFTW
#ifdef NAMD_FFTW_3
#include <fftw3.h>
#else
// fftw2 doesn't have these defined
#define fftwf_malloc fftw_malloc
#define fftwf_free fftw_free
#ifdef NAMD_FFTW_NO_TYPE_PREFIX
#include <fftw.h>
#include <rfftw.h>
#else
#include <sfftw.h>
#include <srfftw.h>
#endif
#endif
#endif
#if defined(WIN32) && !defined(__CYGWIN__)
#include <direct.h>
#define CHDIR _chdir
#define MKDIR(X) mkdir(X)
#define PATHSEP '\\'
#define PATHSEPSTR "\\"
#else
#include <unistd.h>
#define CHDIR chdir
#define MKDIR(X) mkdir(X,0777)
#define PATHSEP '/'
#define PATHSEPSTR "/"
#endif
#include <fcntl.h>
#include <sys/stat.h>
#ifdef WIN32
#include <io.h>
#define access(PATH,MODE) _access(PATH,00)
#endif
#include <fstream>
using namespace std;

#ifdef WIN32
extern "C" {
  double erfc(double);
}
#endif

#include "strlib.h"    //  For strcasecmp and strncasecmp

#include "ComputeNonbondedMICKernel.h"

//#ifdef NAMD_CUDA
//bool one_cuda_device_per_node();
//#endif
#include "DeviceCUDA.h"
#ifdef NAMD_CUDA
#ifdef WIN32
#define __thread __declspec(thread)
#endif
extern __thread DeviceCUDA *deviceCUDA;
#endif

//#define DEBUGM
#include "Debug.h"

#define XXXBIGREAL 1.0e32


char* SimParameters::getfromparseopts(const char* name, char *outbuf) {
  if ( parseopts ) return parseopts->getfromptr(name,outbuf);
  else return 0;
}

int SimParameters::istrueinparseopts(const char* name) {
  if ( parseopts ) return parseopts->istruefromptr(name);
  else return -1;
}

int SimParameters::issetinparseopts(const char* name) {
  if ( parseopts ) return parseopts->issetfromptr(name);
  else return -1;
}


/************************************************************************/
/*                  */
/*      FUNCTION initialize_config_data      */
/*                  */
/*  This function is used by the master process to populate the     */
/*   simulation parameters from a ConfigList object that is passed to   */
/*   it.  Each parameter is checked to make sure that it has a value    */
/*   that makes sense, and that it doesn't conflict with any other      */
/*   values that have been given.          */
/*                  */
/************************************************************************/

void SimParameters::initialize_config_data(ConfigList *config, char *&cwd)

{

   parseopts = new ParseOptions;   //  Object to check consistency of config file
   ParseOptions &opts = *parseopts;

   config_parser(opts);

   ///////////////////////////////// check the internal consistency
   if (!opts.check_consistency())
   {
      NAMD_die("Internal error in configuration file parser");
   }

   // Now, feed the object with the actual configuration options through the
   // ParseOptions file and make sure everything is OK
   if (!opts.set(*config))
   {
      NAMD_die("ERROR(S) IN THE CONFIGURATION FILE");
   }

   //// now do more logic stuff that can't be done by the ParseOptions object

   check_config(opts,config,cwd);

   print_config(opts,config,cwd);

}

/************************************************************************/
/*                                                                      */
/*      FUNCTION scriptSet                                              */
/*                                                                      */
/************************************************************************/

int atobool(const char *s) {
  return ( (! strncasecmp(s,"yes",8)) ||
           (! strncasecmp(s,"on",8)) ||
           (! strncasecmp(s,"true",8)) );
};

void SimParameters::scriptSet(const char *param, const char *value) {

  if ( CkMyRank() ) return;

#define MAX_SCRIPT_PARAM_SIZE 128
#define SCRIPT_PARSE_BOOL(NAME,VAR) { if ( ! strncasecmp(param,(NAME),MAX_SCRIPT_PARAM_SIZE) ) { (VAR) = atobool(value); return; } }
#define SCRIPT_PARSE_INT(NAME,VAR) { if ( ! strncasecmp(param,(NAME),MAX_SCRIPT_PARAM_SIZE) ) { (VAR) = atoi(value); return; } }
#define SCRIPT_PARSE_FLOAT(NAME,VAR) { if ( ! strncasecmp(param,(NAME),MAX_SCRIPT_PARAM_SIZE) ) { (VAR) = atof(value); return; } }
#define SCRIPT_PARSE_MOD_FLOAT(NAME,VAR,MOD) { if ( ! strncasecmp(param,(NAME),MAX_SCRIPT_PARAM_SIZE) ) { (VAR) = atof(value) MOD; return; } }
#define SCRIPT_PARSE_VECTOR(NAME,VAR) { if ( ! strncasecmp(param,(NAME),MAX_SCRIPT_PARAM_SIZE) ) { (VAR).set(value); return; } }
#define SCRIPT_PARSE_STRING(NAME,VAR) { if ( ! strncasecmp(param,(NAME),MAX_SCRIPT_PARAM_SIZE) ) { strcpy(VAR,value); return; } }

  SCRIPT_PARSE_FLOAT("scriptArg1",scriptArg1)
  SCRIPT_PARSE_FLOAT("scriptArg2",scriptArg2)
  SCRIPT_PARSE_FLOAT("scriptArg3",scriptArg3)
  SCRIPT_PARSE_FLOAT("scriptArg4",scriptArg4)
  SCRIPT_PARSE_FLOAT("scriptArg5",scriptArg5)
  SCRIPT_PARSE_INT("scriptIntArg1",scriptIntArg1)
  SCRIPT_PARSE_INT("scriptIntArg2",scriptIntArg2)
  SCRIPT_PARSE_STRING("scriptStringArg1",scriptStringArg1)
  SCRIPT_PARSE_STRING("scriptStringArg2",scriptStringArg2)
  SCRIPT_PARSE_INT("numsteps",N)
  if ( ! strncasecmp(param,"firsttimestep",MAX_SCRIPT_PARAM_SIZE) ) {
    N = firstTimestep = atoi(value); return;
  }
  SCRIPT_PARSE_FLOAT("reassignTemp",reassignTemp)
  SCRIPT_PARSE_FLOAT("rescaleTemp",rescaleTemp)
  SCRIPT_PARSE_BOOL("velocityQuenching",minimizeOn)
  SCRIPT_PARSE_BOOL("maximumMove",maximumMove)
  // SCRIPT_PARSE_BOOL("Langevin",langevinOn)
  if ( ! strncasecmp(param,"Langevin",MAX_SCRIPT_PARAM_SIZE) ) {
    langevinOn = atobool(value);
    if ( langevinOn && ! langevinOnAtStartup ) {
      NAMD_die("Langevin must be enabled at startup to disable and re-enable in script.");
    }
    return;
  }
  SCRIPT_PARSE_FLOAT("langevinTemp",langevinTemp)
  SCRIPT_PARSE_BOOL("langevinBAOAB",langevin_useBAOAB) // [!!] Use the BAOAB integrator or not
  SCRIPT_PARSE_FLOAT("loweAndersenTemp",loweAndersenTemp) // BEGIN LA, END LA
  SCRIPT_PARSE_BOOL("stochRescale",stochRescaleOn)
  SCRIPT_PARSE_FLOAT("stochRescaleTemp",stochRescaleTemp)
  SCRIPT_PARSE_FLOAT("initialTemp",initialTemp)
  SCRIPT_PARSE_BOOL("useGroupPressure",useGroupPressure)
  SCRIPT_PARSE_BOOL("useFlexibleCell",useFlexibleCell)
  SCRIPT_PARSE_BOOL("useConstantArea",useConstantArea)
  SCRIPT_PARSE_BOOL("fixCellDims",fixCellDims)
  SCRIPT_PARSE_BOOL("fixCellDimX",fixCellDimX)
  SCRIPT_PARSE_BOOL("fixCellDimY",fixCellDimY)
  SCRIPT_PARSE_BOOL("fixCellDimZ",fixCellDimZ)
  SCRIPT_PARSE_BOOL("useConstantRatio",useConstantRatio)
  SCRIPT_PARSE_BOOL("LangevinPiston",langevinPistonOn)
  SCRIPT_PARSE_MOD_FLOAT("LangevinPistonTarget",
			langevinPistonTarget,/PRESSUREFACTOR)
  SCRIPT_PARSE_FLOAT("LangevinPistonPeriod",langevinPistonPeriod)
  SCRIPT_PARSE_FLOAT("LangevinPistonDecay",langevinPistonDecay)
  SCRIPT_PARSE_FLOAT("LangevinPistonTemp",langevinPistonTemp)
  SCRIPT_PARSE_MOD_FLOAT("SurfaceTensionTarget",
			surfaceTensionTarget,*(100.0/PRESSUREFACTOR))
  SCRIPT_PARSE_BOOL("BerendsenPressure",berendsenPressureOn)
  SCRIPT_PARSE_MOD_FLOAT("BerendsenPressureTarget",
			berendsenPressureTarget,/PRESSUREFACTOR)
  SCRIPT_PARSE_MOD_FLOAT("BerendsenPressureCompressibility",
			berendsenPressureCompressibility,*PRESSUREFACTOR)
  SCRIPT_PARSE_FLOAT("BerendsenPressureRelaxationTime",
				berendsenPressureRelaxationTime)
  SCRIPT_PARSE_FLOAT("constraintScaling",constraintScaling)
  SCRIPT_PARSE_FLOAT("consForceScaling",consForceScaling)
  SCRIPT_PARSE_BOOL("drudeHardWall",drudeHardWallOn)
  SCRIPT_PARSE_FLOAT("drudeBondConst",drudeBondConst)
  SCRIPT_PARSE_FLOAT("drudeBondLen",drudeBondLen)
  SCRIPT_PARSE_STRING("outputname",outputFilename)
  SCRIPT_PARSE_INT("outputEnergies",outputEnergies)
  SCRIPT_PARSE_STRING("restartname",restartFilename)
  SCRIPT_PARSE_INT("DCDfreq",dcdFrequency)
  if ( ! strncasecmp(param,"DCDfile",MAX_SCRIPT_PARAM_SIZE) ) {
    close_dcdfile();  // *** implemented in Output.C ***
    strcpy(dcdFilename,value);
    return;
  }
  if ( ! strncasecmp(param,"velDCDfile",MAX_SCRIPT_PARAM_SIZE) ) {
    close_veldcdfile();  // *** implemented in Output.C ***
    strcpy(velDcdFilename,value);
    return;
  }
  SCRIPT_PARSE_STRING("tclBCArgs",tclBCArgs)
  SCRIPT_PARSE_VECTOR("eField",eField)
  SCRIPT_PARSE_FLOAT("eFieldFreq",eFieldFreq)
  SCRIPT_PARSE_FLOAT("eFieldPhase",eFieldPhase)
  SCRIPT_PARSE_FLOAT("accelMDE",accelMDE)
  SCRIPT_PARSE_FLOAT("accelMDalpha",accelMDalpha)
  SCRIPT_PARSE_FLOAT("accelMDTE",accelMDTE)
  SCRIPT_PARSE_FLOAT("accelMDTalpha",accelMDTalpha)
  SCRIPT_PARSE_FLOAT("accelMDGSigma0P",accelMDGSigma0P)
  SCRIPT_PARSE_FLOAT("accelMDGSigma0D",accelMDGSigma0D)
  SCRIPT_PARSE_STRING("accelMDGRestartFile",accelMDGRestartFile)
  SCRIPT_PARSE_VECTOR("stirAxis",stirAxis)
  SCRIPT_PARSE_VECTOR("stirPivot",stirPivot)
  if ( ! strncasecmp(param,"mgridforcescale",MAX_SCRIPT_PARAM_SIZE) ) {
    NAMD_die("Can't yet modify mgridforcescale in a script");
    return;
  }
  if ( ! strncasecmp(param,"mgridforcevoff",MAX_SCRIPT_PARAM_SIZE) ) {
    NAMD_die("Can't yet modify mgridforcevoff in a script");
    return;
  }

  if ( ! strncasecmp(param,"fixedatoms",MAX_SCRIPT_PARAM_SIZE) ) {
    if ( ! fixedAtomsOn )
      NAMD_die("FixedAtoms may not be enabled in a script.");
    if ( ! fixedAtomsForces )
      NAMD_die("To use fixedAtoms in script first use fixedAtomsForces yes.");
    fixedAtomsOn = atobool(value);
    return;
  }

//fepb
  if ( ! strncasecmp(param,"alch",MAX_SCRIPT_PARAM_SIZE) ) {
    alchOn = atobool(value);
    if ( alchOn && ! alchOnAtStartup ) {
       NAMD_die("Alchemy must be enabled at startup to disable and re-enable in script.");
    }
    alchFepOn = alchOn && alchFepOnAtStartup;
    alchThermIntOn = alchOn && alchThermIntOnAtStartup;
    ComputeNonbondedUtil::select();
    if ( PMEOn ) ComputePmeUtil::select();
    return;
  }
  SCRIPT_PARSE_INT("alchEquilSteps",alchEquilSteps)

  if ( ! strncasecmp(param,"alchLambda",MAX_SCRIPT_PARAM_SIZE) ) {
    alchLambda = atof(value);
    if ( alchLambda < 0.0 || 1.0 < alchLambda ) {
      NAMD_die("Alchemical lambda values should be in the range [0.0, 1.0]\n");
    }
    ComputeNonbondedUtil::select();
    return;
  }

  if ( ! strncasecmp(param,"alchLambda2",MAX_SCRIPT_PARAM_SIZE) ) {
    alchLambda2 = atof(value);
    if ( alchLambda2 < 0.0 || 1.0 < alchLambda2 ) {
      NAMD_die("Alchemical lambda values should be in the range [0.0, 1.0]\n");
    }
    ComputeNonbondedUtil::select();
    return;
  }

  if ( ! strncasecmp(param,"alchLambdaIDWS",MAX_SCRIPT_PARAM_SIZE) ) {
    alchLambdaIDWS = atof(value);
    setupIDWS();
    ComputeNonbondedUtil::select();
    return;
  }

  if ( ! strncasecmp(param,"alchLambdaFreq",MAX_SCRIPT_PARAM_SIZE) ) {
    alchLambdaFreq = atoi(value);
    if ( alchLambdaIDWS >= 0 ) {
      NAMD_die("alchLambdaIDWS and alchLambdaFreq are not compatible.\n");
    }
    ComputeNonbondedUtil::select();
    return;
  }
//fepe

  if ( ! strncasecmp(param,"nonbondedScaling",MAX_SCRIPT_PARAM_SIZE) ) {
    nonbondedScaling = atof(value);
    ComputeNonbondedUtil::select();
    return;
  }

  if ( ! strncasecmp(param,"commOnly",MAX_SCRIPT_PARAM_SIZE) ) {
    commOnly = atobool(value);
    ComputeNonbondedUtil::select();
    return;
  }

  // REST2 - We don't have to make any changes to ComputeNonbondedUtil
  // other than recalculating the LJTable.  Skip doing the other stuff.
  if ( ! strncasecmp(param,"soluteScalingFactor",MAX_SCRIPT_PARAM_SIZE)) {
    soluteScalingFactor = atof(value);
    if (soluteScalingFactor < 0.0) {
      NAMD_die("Solute scaling factor should be non-negative\n");
    }
    soluteScalingFactorCharge = soluteScalingFactor;
    soluteScalingFactorVdw = soluteScalingFactor;
    // update LJTable for CPU
    ComputeNonbondedUtil::select();
#ifdef NAMD_CUDA
    // update LJTable for GPU, needs CPU update first
    ComputeCUDAMgr::getComputeCUDAMgr()->update();
#endif
    return;
  }
  if ( ! strncasecmp(param,"soluteScalingFactorVdw",MAX_SCRIPT_PARAM_SIZE)) {
    soluteScalingFactorVdw = atof(value);
    if (soluteScalingFactorVdw < 0.0) {
      NAMD_die("Solute scaling factor for van der Waals "
          "should be non-negative\n");
    }
    // update LJTable for CPU
    ComputeNonbondedUtil::select();
#ifdef NAMD_CUDA
    // update LJTable for GPU, needs CPU update first
    ComputeCUDAMgr::getComputeCUDAMgr()->update();
#endif
    return;
  }
  if ( ! strncasecmp(param,"soluteScalingFactorCharge",MAX_SCRIPT_PARAM_SIZE)) {
    soluteScalingFactorCharge = atof(value);
    if (soluteScalingFactorCharge < 0.0) {
      NAMD_die("Solute scaling factor for electrostatics "
          "should be non-negative\n");
    }
    return;
  }

  char *error = new char[2 * MAX_SCRIPT_PARAM_SIZE + 100];
  sprintf(error,"Setting parameter %s from script failed!\n",param);
  NAMD_die(error);

}

void SimParameters::nonbonded_select() {
    ComputeNonbondedUtil::select();
}

void SimParameters::pme_select() {
  ComputePmeUtil::select();
}

/************************************************************************/
/*                                                                      */
/*      FUNCTION config_parser                                          */
/*                                                                      */
/************************************************************************/

void SimParameters::config_parser(ParseOptions &opts) {

   //  Set all variable to fallback default values.  This is not really
   //  necessary, as we give default values when we set up the ParseOptions
   //  object, but it helps the debuggers figure out we've initialized the
   //  variables.
   HydrogenBonds = FALSE;
   useAntecedent = TRUE;
   aaAngleExp = 2;
   haAngleExp = 4;
   distAttExp = 4;
   distRepExp = 6;
   dhaCutoffAngle = 100.0;
   dhaOnAngle = 60.0;
   dhaOffAngle = 80.0;
   daCutoffDist = 7.5;
   daOnDist = 5.5;
   daOffDist = 6.5;

   config_parser_basic(opts);
   config_parser_fileio(opts);
   config_parser_fullelect(opts);
   config_parser_methods(opts);
   config_parser_constraints(opts);
   #ifdef OPENATOM_VERSION
   config_parser_openatom(opts);
   #endif // OPENATOM_VERSION

   config_parser_gridforce(opts);
   config_parser_mgridforce(opts);
   config_parser_movdrag(opts);
   config_parser_rotdrag(opts);
   config_parser_constorque(opts);
   config_parser_boundary(opts);
   config_parser_misc(opts);

}

void SimParameters::config_parser_basic(ParseOptions &opts) {

   //  So first we set up the ParseOptions objects so that it will check
   //  all of the logical rules that the configuration file must follow.

   opts.optional("main", "obsolete", "used to flag obsolete options",
    PARSE_STRING);

   ////// basic options
   opts.require("main", "timestep", "size of the timestep, in fs",
    &dt, 1.0);
   opts.range("timestep", NOT_NEGATIVE);
   opts.units("timestep", N_FSEC);

   opts.optional("main", "numsteps", "number of timesteps to perform",
    &N,0);
   opts.range("numsteps", NOT_NEGATIVE);

   opts.optional("main", "stepspercycle",
      "Number of steps between atom migrations",
      &stepsPerCycle, 20);
   opts.range("stepspercycle", POSITIVE);

   opts.require("main", "cutoff", "local electrostatic and Vdw distance",
      &cutoff);
   opts.range("cutoff", POSITIVE);
   opts.units("cutoff", N_ANGSTROM);

   opts.optional("main", "nonbondedScaling", "nonbonded scaling factor",
     &nonbondedScaling, 1.0);
   opts.range("nonbondedScaling", NOT_NEGATIVE);

   opts.optional("main", "limitDist", "limit nonbonded below this distance",
     &limitDist, 0.0);
   opts.range("limitDist", NOT_NEGATIVE);

   opts.require("main", "exclude", "Electrostatic and VDW exclusion policy",
    PARSE_STRING);

   opts.optional("exclude", "1-4scaling", "1-4 electrostatic scaling factor",
     &scale14, 1.0);
   opts.range("1-4scaling", POSITIVE);

   opts.optionalB("main", "switching",
     "Should a smoothing function be used?", &switchingActive, TRUE);

   opts.optionalB("switching", "vdwForceSwitching",
     "Use force switching for vdw?", &vdwForceSwitching, FALSE);

   opts.optional("switching", "switchdist",
     "Distance for switching function activation",
     &switchingDist);
   opts.range("switchdist", POSITIVE);
   opts.units("switchdist", N_ANGSTROM);

   opts.optionalB("main", "martiniSwitching",
     "Use Martini residue-based coarse-grain switching?", &martiniSwitching, FALSE);
   opts.optionalB("main", "martiniDielAllow",
     "Allow use of dielectric != 15.0 when using Martini", &martiniDielAllow, FALSE);

   opts.optional("main", "pairlistdist",  "Pairlist inclusion distance",
     &pairlistDist);
   opts.range("pairlistdist", POSITIVE);
   opts.units("pairlistdist", N_ANGSTROM);

   opts.optional("main", "pairlistMinProcs",  "Min procs for pairlists",
     &pairlistMinProcs,1);
   opts.range("pairlistMinProcs", POSITIVE);

   opts.optional("main", "pairlistsPerCycle",  "regenerate x times per cycle",
     &pairlistsPerCycle,2);
   opts.range("pairlistsPerCycle", POSITIVE);

   opts.optional("main", "outputPairlists", "how often to print warnings",
     &outputPairlists, 0);
   opts.range("outputPairlists", NOT_NEGATIVE);

   opts.optional("main", "pairlistShrink",  "tol *= (1 - x) on regeneration",
     &pairlistShrink,0.01);
   opts.range("pairlistShrink", NOT_NEGATIVE);

   opts.optional("main", "pairlistGrow",  "tol *= (1 + x) on trigger",
     &pairlistGrow, 0.01);
   opts.range("pairlistGrow", NOT_NEGATIVE);

   opts.optional("main", "pairlistTrigger",  "trigger is atom > (1 - x) * tol",
     &pairlistTrigger, 0.3);
   opts.range("pairlistTrigger", NOT_NEGATIVE);

   opts.optional("main", "temperature", "initial temperature",
     &initialTemp);
   opts.range("temperature", NOT_NEGATIVE);
   opts.units("temperature", N_KELVIN);

   opts.optionalB("main", "COMmotion", "allow initial center of mass movement",
      &comMove, FALSE);

   opts.optionalB("main", "zeroMomentum", "constrain center of mass",
      &zeroMomentum, FALSE);
   opts.optionalB("zeroMomentum", "zeroMomentumAlt", "constrain center of mass",
      &zeroMomentumAlt, FALSE);

   opts.optionalB("main", "wrapWater", "wrap waters around periodic boundaries on output",
      &wrapWater, FALSE);
   opts.optionalB("main", "wrapAll", "wrap all clusters around periodic boundaries on output",
      &wrapAll, FALSE);
   opts.optionalB("main", "wrapNearest", "wrap to nearest image to cell origin",
      &wrapNearest, FALSE);

   opts.optional("main", "dielectric", "dielectric constant",
     &dielectric, 1.0);
   opts.range("dielectric", POSITIVE); // Hmmm, dielectric < 1 ...

   opts.optional("main", "margin", "Patch width margin", &margin, XXXBIGREAL);
   opts.range("margin", NOT_NEGATIVE);
   opts.units("margin", N_ANGSTROM);

   opts.optional("main", "seed", "Initial random number seed", &randomSeed);
   opts.range("seed", POSITIVE);

   opts.optional("main", "outputEnergies", "How often to print energies in timesteps",
     &outputEnergies, 1);
   opts.range("outputEnergies", POSITIVE);

   opts.optional("main", "outputMomenta", "How often to print linear and angular momenta in timesteps",
     &outputMomenta, 0);
   opts.range("outputMomenta", NOT_NEGATIVE);

   opts.optional("main", "outputTiming", "How often to print timing data in timesteps",
     &outputTiming);
   opts.range("outputTiming", NOT_NEGATIVE);

   opts.optional("main", "outputCudaTiming", "How often to print CUDA timing data in timesteps",
     &outputCudaTiming, 0);
   opts.range("outputCudaTiming", NOT_NEGATIVE);

   opts.optional("main", "outputPressure", "How often to print pressure data in timesteps",
     &outputPressure, 0);
   opts.range("outputPressure", NOT_NEGATIVE);

   opts.optionalB("main", "mergeCrossterms", "merge crossterm energy with dihedral when printing?",
      &mergeCrossterms, TRUE);

   opts.optional("main", "MTSAlgorithm", "Multiple timestep algorithm",
    PARSE_STRING);

   opts.optional("main", "longSplitting", "Long range force splitting option",
    PARSE_STRING);

   opts.optionalB("main", "ignoreMass", "Do not use masses to find hydrogen atoms",
    &ignoreMass, FALSE);

   opts.optional("main", "splitPatch", "Atom into patch splitting option",
    PARSE_STRING);
   opts.optional("main", "hgroupCutoff", "Hydrogen margin", &hgroupCutoff, 2.5);

   opts.optional("main", "extendedSystem",
    "Initial configuration of extended system variables and periodic cell",
    PARSE_STRING);

   opts.optional("main", "cellBasisVector1", "Basis vector for periodic cell",
    &cellBasisVector1);
   opts.optional("main", "cellBasisVector2", "Basis vector for periodic cell",
    &cellBasisVector2);
   opts.optional("main", "cellBasisVector3", "Basis vector for periodic cell",
    &cellBasisVector3);
   opts.optional("main", "cellOrigin", "Fixed center of periodic cell",
    &cellOrigin);

   opts.optionalB("main", "molly", "Rigid bonds to hydrogen",&mollyOn,FALSE);
   opts.optional("main", "mollyTolerance", "Error tolerance for MOLLY",
                 &mollyTol, 0.00001);
   opts.optional("main", "mollyIterations",
		 "Max number of iterations for MOLLY", &mollyIter, 100);

   opts.optional("main", "rigidBonds", "Rigid bonds to hydrogen",PARSE_STRING);
   opts.optional("main", "rigidTolerance",
                 "Error tolerance for rigid bonds to hydrogen",
                 &rigidTol, 1.0e-8);
   opts.optional("main", "rigidIterations",
		 "Max number of SHAKE iterations for rigid bonds to hydrogen",
		 &rigidIter, 100);
   opts.optionalB("main", "rigidDieOnError",
		 "Die if rigidTolerance is not achieved after rigidIterations",
		 &rigidDie, TRUE);
   opts.optionalB("main", "useSettle",
                  "Use the SETTLE algorithm for rigid waters",
                 &useSettle, TRUE);

   opts.optional("main", "nonbondedFreq", "Nonbonded evaluation frequency",
    &nonbondedFrequency, 1);
   opts.range("nonbondedFreq", POSITIVE);

   opts.optionalB("main", "outputPatchDetails", "print number of atoms in each patch",
      &outputPatchDetails, FALSE);
   opts.optionalB("main", "staticAtomAssignment", "never migrate atoms",
      &staticAtomAssignment, FALSE);
   opts.optionalB("main", "replicaUniformPatchGrids", "same patch grid size on all replicas",
      &replicaUniformPatchGrids, FALSE);
#ifndef MEM_OPT_VERSION
   // in standard (non-mem-opt) version, enable lone pairs by default
   // for compatibility with recent force fields
   opts.optionalB("main", "lonePairs", "Enable lone pairs", &lonepairs, TRUE);
#else
   // in mem-opt version, disable lone pairs by default
   // because they are not supported
   opts.optionalB("main", "lonePairs", "Enable lone pairs", &lonepairs, FALSE);
#endif
   opts.optional("main", "waterModel", "Water model to use", PARSE_STRING);
   opts.optionalB("main", "LJcorrection", "Apply analytical tail corrections for energy and virial", &LJcorrection, FALSE);
#ifdef TIMER_COLLECTION
   opts.optional("main", "TimerBinWidth",
       "Bin width of timer histogram collection in microseconds",
       &timerBinWidth, 1.0);
#endif
#if defined(NAMD_NVTX_ENABLED) || defined(NAMD_CMK_TRACE_ENABLED)
   // default NVTX or Projections profiling is up to the first 1000 patches
   opts.optional("main", "beginEventPatchID","Beginning patch ID for profiling",
       &beginEventPatchID, 0);
   opts.optional("main", "endEventPatchID", "Ending patch ID for profiling",
       &endEventPatchID, 5000);
   // default NVTX or Projections profiling is up to the first 1000 time steps
   opts.optional("main", "beginEventStep", "Beginning time step for profiling",
       &beginEventStep, 0);
   opts.optional("main", "endEventStep", "Ending time step for profiling",
       &endEventStep, 1000);
#endif
}

void SimParameters::config_parser_fileio(ParseOptions &opts) {

   /////////////// file I/O

   opts.optional("main", "cwd", "current working directory", PARSE_STRING);

// In order to include AMBER options, "coordinates", "structure"
// and "parameters" are now optional, not required. The presence
// of them will be checked later in check_config()

//   opts.require("main", "coordinates", "initial PDB coordinate file",
//    PARSE_STRING);
   opts.optional("main", "coordinates", "initial PDB coordinate file",
    PARSE_STRING);

   opts.optional("main", "velocities",
     "initial velocities, given as a PDB file", PARSE_STRING);
   opts.optional("main", "binvelocities",
     "initial velocities, given as a binary restart", PARSE_STRING);
   opts.optional("main", "bincoordinates",
     "initial coordinates in a binary restart file", PARSE_STRING);
#ifdef MEM_OPT_VERSION
   opts.optional("main", "binrefcoords",
     "reference coordinates in a binary restart file", PARSE_STRING);
#endif

//   opts.require("main", "structure", "initial PSF structure file",
//    PARSE_STRING);
   opts.optional("main", "structure", "initial PSF structure file",
    PARSE_STRING);

//   opts.require("main", "parameters",
//"CHARMm 19 or CHARMm 22 compatable force field file (multiple "
//"inputs allowed)", PARSE_MULTIPLES);
   opts.optional("main", "parameters",
"CHARMm 19 or CHARMm 22 compatable force field file (multiple "
"inputs allowed)", PARSE_MULTIPLES);


   //****** BEGIN CHARMM/XPLOR type changes
   //// enable XPLOR as well as CHARMM input files for parameters
   opts.optionalB("parameters", "paraTypeXplor", "Parameter file in Xplor format?", &paraTypeXplorOn, FALSE);
   opts.optionalB("parameters", "paraTypeCharmm", "Parameter file in Charmm format?", &paraTypeCharmmOn, FALSE);
   //****** END CHARMM/XPLOR type changes

   // Ported by JLai -- JE - Go parameters
   opts.optionalB("main", "GromacsPair", "Separately calculate pair interactions", &goGroPair, FALSE);
   opts.optionalB("main", "GoForcesOn", "Go forces will be calculated", &goForcesOn, FALSE);
   opts.require("GoForcesOn", "GoParameters", "Go parameter file", goParameters);
   opts.require("GoForcesOn", "GoCoordinates", "target coordinates for Go forces", goCoordinates);
   // Added by JLai -- Go-method parameter -- switches between using the matrix and sparse matrix representations [6.3.11]
   //   opts.optional("GoForcesOn", "GoMethod", "Which type of matrix should be used to store Go contacts?", &goMethod);
   // Added by JLai -- Go-method parameter [8.2.11]
   //opts.optional("GoForcesOn", "GoMethod", "Which type of matrix should be used to store Go contacts?", PARSE_STRING);
   opts.require("GoForcesOn", "GoMethod", "Which type of matrix should be used to store Go contacts?", PARSE_STRING);
  // End of Port -- JL

   opts.require("main", "outputname",
    "prefix for the final PDB position and velocity filenames",
    outputFilename);

   opts.optional("main", "auxFile", "Filename for data stream output",
     auxFilename);

   opts.optional("main", "numinputprocs", "Number of pes to use for parallel input",
                 &numinputprocs, 0);
   opts.range("numinputprocs", NOT_NEGATIVE);

   opts.optional("main", "numoutputprocs", "Number of pes to use for parallel output",
                 &numoutputprocs, 0);
   opts.range("numoutputprocs", NOT_NEGATIVE);
   opts.optional("main", "numoutputwriters", "Number of output processors that simultaneously write to an output file",
                 &numoutputwrts, 1);
   opts.range("numoutputwriters", NOT_NEGATIVE);

   opts.optional("main", "DCDfreq", "Frequency of DCD trajectory output, in "
    "timesteps", &dcdFrequency, 0);
   opts.range("DCDfreq", NOT_NEGATIVE);
   opts.optional("DCDfreq", "DCDfile", "DCD trajectory output file name",
     dcdFilename);
   opts.optionalB("DCDfreq", "DCDunitcell", "Store unit cell in dcd timesteps?",
       &dcdUnitCell);

   opts.optional("main", "velDCDfreq", "Frequency of velocity "
    "DCD output, in timesteps", &velDcdFrequency, 0);
   opts.range("velDCDfreq", NOT_NEGATIVE);
   opts.optional("velDCDfreq", "velDCDfile", "velocity DCD output file name",
     velDcdFilename);

   opts.optional("main", "forceDCDfreq", "Frequency of force"
    "DCD output, in timesteps", &forceDcdFrequency, 0);
   opts.range("forceDCDfreq", NOT_NEGATIVE);
   opts.optional("forceDCDfreq", "forceDCDfile", "force DCD output file name",
     forceDcdFilename);

   opts.optional("main", "XSTfreq", "Frequency of XST trajectory output, in "
    "timesteps", &xstFrequency, 0);
   opts.range("XSTfreq", NOT_NEGATIVE);
   opts.optional("XSTfreq", "XSTfile", "Extended sytem trajectory output "
    "file name", xstFilename);

   opts.optional("main", "restartfreq", "Frequency of restart file "
    "generation", &restartFrequency, 0);
   opts.range("restartfreq", NOT_NEGATIVE);
   opts.optional("restartfreq", "restartname", "Prefix for the position and "
     "velocity PDB files used for restarting", restartFilename);
   opts.optionalB("restartfreq", "restartsave", "Save restart files with "
     "unique filenames rather than overwriting", &restartSave, FALSE);
   opts.optionalB("restartfreq", "restartsavedcd", "Save DCD files with "
     "unique filenames at each restart", &restartSaveDcd, FALSE);

   opts.optionalB("restartfreq", "binaryrestart", "Specify use of binary restart files ",
       &binaryRestart, TRUE);

   opts.optionalB("outputname", "binaryoutput", "Specify use of binary output files ",
       &binaryOutput, TRUE);

   opts.optionalB("main", "amber", "Is it AMBER force field?",
       &amberOn, FALSE);
   opts.optionalB("amber", "readexclusions", "Read exclusions from parm file?",
       &readExclusions, TRUE);
   opts.require("amber", "scnb", "1-4 VDW interactions are divided by scnb",
       &vdwscale14, 2.0);
   opts.require("amber", "parmfile", "AMBER parm file", PARSE_STRING);
   opts.optional("amber", "ambercoor", "AMBER coordinate file", PARSE_STRING);

   /* GROMACS options */
   opts.optionalB("main", "gromacs", "Use GROMACS-like force field?",
       &gromacsOn, FALSE);
   opts.require("gromacs", "grotopfile", "GROMACS topology file",
		PARSE_STRING);
   opts.optional("gromacs", "grocoorfile","GROMACS coordinate file",
		 PARSE_STRING);

  // OPLS options
   opts.optionalB("main", "vdwGeometricSigma",
       "Use geometric mean to combine L-J sigmas, as for OPLS",
       &vdwGeometricSigma, FALSE);

   // load/store computeMap
   opts.optional("main", "computeMapFile", "Filename for computeMap",
     computeMapFilename);
   opts.optionalB("main", "storeComputeMap", "store computeMap?",
       &storeComputeMap, FALSE);
   opts.optionalB("main", "loadComputeMap", "load computeMap?",
       &loadComputeMap, FALSE);
}


void SimParameters::config_parser_fullelect(ParseOptions &opts) {

   /////////// FMA options
#ifdef DPMTA
   DebugM(1,"DPMTA setup start\n");
   //  PMTA is included, so really get these values
   opts.optionalB("main", "FMA", "Should FMA be used?", &FMAOn, FALSE);
   opts.optional("FMA", "FMALevels", "Tree levels to use in FMA", &FMALevels,
     5);
   opts.range("FMALevels", POSITIVE);
   opts.optional("FMA", "FMAMp", "Number of FMA multipoles", &FMAMp, 8);
   opts.range("FMAMp", POSITIVE);
   opts.optionalB("FMA", "FMAFFT", "Use FFT enhancement in FMA?", &FMAFFTOn, TRUE);
   opts.optional("FMAFFT", "FMAFFTBlock", "FFT blocking factor",
    &FMAFFTBlock, 4);
   opts.range("FMAFFTBlock", POSITIVE);
   DebugM(1,"DPMTA setup end\n");
#else
   //  PMTA is NOT included.  So just set all the values to 0.
   FMAOn = FALSE;
   FMALevels = 0;
   FMAMp = 0;
   FMAFFTOn = FALSE;
   FMAFFTBlock = 0;
#endif

   opts.optional("main", "fullElectFrequency",
      "Number of steps between full electrostatic executions",
      &fullElectFrequency);
   opts.range("fullElectFrequency", POSITIVE);

   //  USE OF THIS PARAMETER DISCOURAGED
   opts.optional("main", "fmaFrequency",
      "Number of steps between full electrostatic executions",
      &fmaFrequency);
   opts.range("fmaFrequency", POSITIVE);

   opts.optional("main", "fmaTheta",
      "FMA theta parameter value",
      &fmaTheta,0.715);
   opts.range("fmaTheta", POSITIVE);

   opts.optionalB("main", "FullDirect", "Should direct calculations of full electrostatics be performed?",
      &fullDirectOn, FALSE);


   ///////////  Multilevel Summation Method

   opts.optionalB("main", "MSM",
       "Use multilevel summation method for electrostatics?",
       &MSMOn, FALSE);
   opts.optional("MSM", "MSMQuality", "MSM quality",
       &MSMQuality, 0);
   opts.optional("MSM", "MSMApprox", "MSM approximation",
       &MSMApprox, 0);
   opts.optional("MSM", "MSMSplit", "MSM splitting",
       &MSMSplit, 0);
   opts.optional("MSM", "MSMLevels", "MSM maximum number of levels",
       &MSMLevels, 0);  // set to 0 adapts to as many as needed
   opts.optional("MSM", "MSMGridSpacing", "MSM grid spacing (Angstroms)",
       &MSMGridSpacing, 2.5);
   opts.optional("MSM", "MSMPadding", "MSM padding (Angstroms)",
       &MSMPadding, 2.5);
   opts.optional("MSM", "MSMxmin", "MSM x minimum (Angstroms)", &MSMxmin, 0);
   opts.optional("MSM", "MSMxmax", "MSM x maximum (Angstroms)", &MSMxmax, 0);
   opts.optional("MSM", "MSMymin", "MSM y minimum (Angstroms)", &MSMymin, 0);
   opts.optional("MSM", "MSMymax", "MSM y maximum (Angstroms)", &MSMymax, 0);
   opts.optional("MSM", "MSMzmin", "MSM z minimum (Angstroms)", &MSMzmin, 0);
   opts.optional("MSM", "MSMzmax", "MSM z maximum (Angstroms)", &MSMzmax, 0);
   opts.optional("MSM", "MSMBlockSizeX",
       "MSM grid block size along X direction (for decomposing parallel work)",
       &MSMBlockSizeX, 8);
   opts.optional("MSM", "MSMBlockSizeY",
       "MSM grid block size along Y direction (for decomposing parallel work)",
       &MSMBlockSizeY, 8);
   opts.optional("MSM", "MSMBlockSizeZ",
       "MSM grid block size along Z direction (for decomposing parallel work)",
       &MSMBlockSizeZ, 8);

   opts.optionalB("MSM", "MsmSerial",
       "Use MSM serial version for long-range calculation?",
       &MsmSerialOn, FALSE);


   ///////////  Fast Multipole Method

   opts.optionalB("main", "FMM",
       "Use fast multipole method for electrostatics?",
       &FMMOn, FALSE);
   opts.optional("FMM", "FMMLevels", "FMM number of levels",
       &FMMLevels, 0);
   opts.optional("FMM", "FMMPadding", "FMM padding margin (Angstroms)",
       &FMMPadding, 0);

   opts.optionalB("main", "useCUDA2", "Use new CUDA code", &useCUDA2, TRUE);

   ///////////  Particle Mesh Ewald

   opts.optionalB("main", "PME", "Use particle mesh Ewald for electrostatics?",
	&PMEOn, FALSE);
   opts.optional("PME", "PMETolerance", "PME direct space tolerance",
	&PMETolerance, 1.e-6);
   opts.optional("PME", "PMEInterpOrder", "PME interpolation order",
	&PMEInterpOrder, 4);  // cubic interpolation is default
   opts.optional("PME", "PMEGridSizeX", "PME grid in x dimension",
	&PMEGridSizeX, 0);
   opts.optional("PME", "PMEGridSizeY", "PME grid in y dimension",
	&PMEGridSizeY, 0);
   opts.optional("PME", "PMEGridSizeZ", "PME grid in z dimension",
	&PMEGridSizeZ, 0);
   opts.optional("PME", "PMEGridSpacing", "Maximum PME grid spacing (Angstroms)",
	&PMEGridSpacing, 0.);
   opts.range("PMEGridSpacing", NOT_NEGATIVE);
   opts.optional("PME", "PMEProcessors",
	"PME FFT and reciprocal sum processor count", &PMEProcessors, 0);
   opts.optional("PME", "PMEMinSlices",
	"minimum thickness of PME reciprocal sum slab", &PMEMinSlices, 2);
   opts.range("PMEMinSlices", NOT_NEGATIVE);
   opts.optional("PME", "PMEPencils",
	"PME FFT and reciprocal sum pencil grid size", &PMEPencils, -1);
   opts.optional("PME", "PMEPencilsX",
	"PME FFT and reciprocal sum pencil grid size X", &PMEPencilsX, 0);
   opts.optional("PME", "PMEPencilsY",
	"PME FFT and reciprocal sum pencil grid size Y", &PMEPencilsY, 0);
   opts.optional("PME", "PMEPencilsZ",
	"PME FFT and reciprocal sum pencil grid size Z", &PMEPencilsZ, 0);
   opts.range("PMEPencilsX", NOT_NEGATIVE);
   opts.range("PMEPencilsY", NOT_NEGATIVE);
   opts.range("PMEPencilsZ", NOT_NEGATIVE);
   opts.optional("PME", "PMEPencilsYLayout",
	"PME FFT and reciprocal sum Y pencil layout strategy", &PMEPencilsYLayout, 0);
   opts.optional("PME", "PMEPencilsXLayout",
	"PME FFT and reciprocal sum X pencil layout strategy", &PMEPencilsXLayout, 1);
   opts.range("PMEPencilsYLayout", NOT_NEGATIVE);
   opts.range("PMEPencilsXLayout", NOT_NEGATIVE);
   opts.optional("PME", "PMESendOrder",
	"PME message ordering control", &PMESendOrder, 0);
   opts.range("PMESendOrder", NOT_NEGATIVE);
   opts.optional("PME", "PMEMinPoints",
	"minimum points per PME reciprocal sum pencil", &PMEMinPoints, 10000);
   opts.range("PMEMinPoints", NOT_NEGATIVE);
   opts.optionalB("main", "PMEBarrier", "Use barrier in PME?",
	&PMEBarrier, FALSE);
   opts.optionalB("main", "PMEOffload", "Offload PME to accelerator?",
	&PMEOffload);

   opts.optionalB("PME", "usePMECUDA", "Use the PME CUDA version", &usePMECUDA, CmiNumPhysicalNodes() < 5);

#ifdef DPME
   opts.optionalB("PME", "useDPME", "Use old DPME code?", &useDPME, FALSE);
#else
   useDPME = 0;
#endif
   opts.optionalB("main", "FFTWPatient", "Use intensive plan creation to optimize FFTW?",
#ifdef WIN32
	&FFTWPatient, TRUE);
#else
	&FFTWPatient, FALSE);
#endif

   opts.optionalB("main", "FFTWEstimate", "Use estimates to optimize FFTW?",
#ifdef WIN32
	&FFTWEstimate, TRUE);
#else
	&FFTWEstimate, FALSE);
#endif
   opts.optionalB("main", "FFTWUseWisdom", "Read/save wisdom file for FFTW?",
#ifdef WIN32
	&FFTWUseWisdom, FALSE);
#else
	&FFTWUseWisdom, TRUE);
#endif
   opts.optional("FFTWUseWisdom", "FFTWWisdomFile", "File for FFTW wisdom",
	FFTWWisdomFile);

}

void SimParameters::config_parser_methods(ParseOptions &opts) {

   /////////// Special Dynamics Methods
   opts.optionalB("main", "minimization", "Should minimization be performed?",
      &minimizeCGOn, FALSE);
   opts.optionalB("main", "minVerbose", "Print extra minimization diagnostics?",
      &minVerbose, FALSE);
   opts.optional("main", "minTinyStep", "very first minimization steps",
      &minTinyStep, 1.0e-6);
   opts.range("minTinyStep", POSITIVE);
   opts.optional("main", "minBabyStep", "initial minimization steps",
      &minBabyStep, 1.0e-2);
   opts.range("minBabyStep", POSITIVE);
   opts.optional("main", "minLineGoal", "line minimization gradient reduction",
      &minLineGoal, 1.0e-3);
   opts.range("minLineGoal", POSITIVE);

   opts.optionalB("main", "velocityQuenching",
      "Should old-style minimization be performed?", &minimizeOn, FALSE);

   opts.optional("main", "maximumMove", "Maximum atom movement per step", &maximumMove, 0.0);
   opts.range("maximumMove", NOT_NEGATIVE);
   opts.units("maximumMove", N_ANGSTROM);

   opts.optionalB("main", "Langevin", "Should Langevin dynamics be performed?",
      &langevinOn, FALSE);
   opts.require("Langevin", "langevinTemp", "Temperature for heat bath in Langevin "
     "dynamics", &langevinTemp);
   opts.range("langevinTemp", NOT_NEGATIVE);
   opts.units("langevinTemp", N_KELVIN);
   opts.optional("Langevin", "langevinDamping", "Damping coefficient (1/ps)",
      &langevinDamping);
   opts.range("langevinDamping", POSITIVE);
   opts.optionalB("Langevin", "langevinHydrogen", "Should Langevin dynamics be applied to hydrogen atoms?",
      &langevinHydrogen);
   opts.optional("Langevin", "langevinFile", "PDB file with temperature "
     "coupling terms (B(i)) (default is the PDB input file)",
     PARSE_STRING);
   opts.optional("Langevin", "langevinCol", "Column in the langevinFile "
     "containing the temperature coupling term B(i);\n"
     "default is 'O'", PARSE_STRING);

   // use BAOAB integration instead of BBK
   opts.optionalB("Langevin", "langevinBAOAB",
       "Should Langevin dynamics be performed using BAOAB integration?",
       &langevin_useBAOAB, FALSE);

// BEGIN LA
   opts.optionalB("main", "LoweAndersen", "Should Lowe-Andersen dynamics be performed?",
		  &loweAndersenOn, FALSE);
   opts.require("LoweAndersen", "loweAndersenTemp", "Temperature for heat bath in Lowe-Andersen "
		"dynamics", &loweAndersenTemp);
   opts.range("loweAndersenTemp", NOT_NEGATIVE);
   opts.units("loweAndersenTemp", N_KELVIN);
   opts.optional("LoweAndersen", "loweAndersenRate", "Collision rate (1/ps)",
		 &loweAndersenRate, 50);
   opts.range("loweAndersenRate", POSITIVE);
   opts.optional("LoweAndersen", "loweAndersenCutoff", "Cutoff radius",
		 &loweAndersenCutoff, 2.7);
   opts.range("loweAndersenCutoff", POSITIVE);
   opts.units("loweAndersenCutoff", N_ANGSTROM);
// END LA

//fepb
   opts.optionalB("main", "alch", "Is achemical simulation being performed?",
     &alchOn, FALSE);
   opts.require("alch", "alchLambda", "Alchemical coupling parameter value",
     &alchLambda);
   opts.range("alchLambda", NOT_NEGATIVE);

   opts.optionalB("alch", "singleTopology",
     "Is single topology used for relative free energy?", &singleTopology, FALSE);

   opts.optionalB("alch", "sdBondScaling",
     "Is S-D bonded terms scaling for relative free energy?", &sdScaling, FALSE);
   
   opts.optional("alch", "alchFile", "PDB file with perturbation flags "
     "default is the input PDB file", PARSE_STRING);
   opts.optional("alch", "alchCol", "Column in the alchFile with the "
     "perturbation flag", PARSE_STRING);

   opts.optional("alch", "unperturbedBondFile", "mini psf file with unperturbed bond info"
     " ", PARSE_STRING);

   opts.optional("alch", "alchOutFreq", "Frequency of alchemical energy"
     "output in timesteps", &alchOutFreq, 5);
   opts.range("alchoutfreq", NOT_NEGATIVE);
   opts.optional("alch", "alchOutFile", "Alchemical energy output filename",
     alchOutFile);

   // soft-core parameters
   opts.optional("alch", "alchVdwShiftCoeff", "Coeff used for generating"
     "the altered alchemical vDW interactions", &alchVdwShiftCoeff, 5.);
   opts.range("alchVdwShiftCoeff", NOT_NEGATIVE);

   opts.optionalB("alch", "alchWCA", "Is WCA decomposition being performed?",
     &alchWCAOn, FALSE);

   // scheduling options for different interaction types
   opts.optional("alch", "alchElecLambdaStart", "Lambda at which electrostatic"
      "scaling of exnihilated particles begins", &alchElecLambdaStart, 0.5);
   opts.range("alchElecLambdaStart", NOT_NEGATIVE);

   opts.optional("alch", "alchVdwLambdaEnd", "Lambda at which vdW"
      "scaling of exnihilated particles ends", &alchVdwLambdaEnd, 1.0);
   opts.range("alchVdwLambdaEnd", NOT_NEGATIVE);

   opts.optional("alch", "alchRepLambdaEnd", "Lambda at which repulsive vdW"
      "scaling of exnihilated particles ends and attractive vdW scaling"
      "begins", &alchRepLambdaEnd, 0.5);
   opts.range("alchRepLambdaEnd", NOT_NEGATIVE);

   opts.optional("alch", "alchBondLambdaEnd", "Lambda at which bonded"
      "scaling of exnihilated particles begins", &alchBondLambdaEnd, 0.0);
   opts.range("alchBondLambdaEnd", NOT_NEGATIVE);

   opts.optionalB("alch", "alchDecouple", "Enable alchemical decoupling?",
     &alchDecouple, FALSE);
   opts.optionalB("alch", "alchBondDecouple", "Enable decoupling of purely "
     "alchemical bonds?", &alchBondDecouple, FALSE);

   // parameters for alchemical analysis options
   opts.optional("alch", "alchType", "Which alchemical method to use?",
     PARSE_STRING);
   opts.optional("alch", "alchLambda2", "Alchemical coupling comparison value",
     &alchLambda2, -1);
   opts.optional("alch", "alchLambdaIDWS", "Alchemical coupling comparison value for interleaved double-wide sampling",
     &alchLambdaIDWS, -1);
   opts.optional("alch", "alchLambdaFreq",
     "Frequency of increasing coupling parameter value", &alchLambdaFreq, 0);
   opts.range("alchLambdaFreq", NOT_NEGATIVE);
   opts.optional("alch", "alchSwitchType", "Switching type flag",
     PARSE_STRING);
   opts.optional("alch", "alchEquilSteps", "Equilibration steps, before "
     "data collection in the alchemical window", &alchEquilSteps, 0);
   opts.range("alchEquilSteps", NOT_NEGATIVE);

   opts.optionalB("alch", "alchEnsembleAvg", "Ensemble Average in use?",
     &alchEnsembleAvg, TRUE);
//fepe

   opts.optionalB("main", "les", "Is locally enhanced sampling enabled?",
     &lesOn, FALSE);
   opts.require("les", "lesFactor", "Local enhancement factor", &lesFactor);
   opts.optional("les", "lesFile", "PDB file with enhancement flags "
     "default is the input PDB file", PARSE_STRING);
   opts.optional("les", "lesCol", "Column in the lesFile with the "
     "enhancement flag", PARSE_STRING);
   opts.optionalB("les", "lesReduceTemp", "Reduce enhanced atom temperature?",
     &lesReduceTemp, FALSE);
   opts.optionalB("les", "lesReduceMass", "Reduce enhanced atom mass?",
     &lesReduceMass, FALSE);

   // REST2 (replica exchange solute tempering) parameters
   opts.optionalB("main", "soluteScaling",
       "Is replica exchange solute tempering enabled?",
       &soluteScalingOn, FALSE);
   opts.require("soluteScaling", "soluteScalingFactor",
       "Solute scaling factor",
       &soluteScalingFactor);
   opts.range("soluteScalingFactor", NOT_NEGATIVE);
   opts.optional("soluteScaling", "soluteScalingFactorCharge",
       "Solute scaling factor for electrostatic interactions",
       &soluteScalingFactorCharge);
   opts.range("soluteScalingFactorCharge", NOT_NEGATIVE);
   opts.optional("soluteScaling", "soluteScalingFactorVdw",
       "Solute scaling factor for van der Waals interactions",
       &soluteScalingFactorVdw);
   opts.range("soluteScalingFactorVdw", NOT_NEGATIVE);
   opts.optional("soluteScaling", "soluteScalingFile",
       "PDB file with scaling flags; if undefined, defaults to main PDB file",
       PARSE_STRING);
   opts.optional("soluteScaling", "soluteScalingCol",
       "Column in the soluteScalingFile providing the scaling flag",
       PARSE_STRING);
   opts.optionalB("main", "soluteScalingAll",
       "Apply scaling also to bond and angle interactions?",
       &soluteScalingAll, FALSE);

   // Drude oscillators
   opts.optionalB("main", "drude", "Perform integration of Drude oscillators?",
       &drudeOn, FALSE);
   opts.require("drude", "drudeTemp", "Temperature for freezing "
       "Drude oscillators", &drudeTemp);
   opts.range("drudeTemp", NOT_NEGATIVE);
   opts.units("drudeTemp", N_KELVIN);
   opts.optional("drude", "drudeDamping", "Damping coefficient (1/ps) for "
       "Drude oscillators", &drudeDamping);
   opts.range("drudeDamping", POSITIVE);
   opts.optional("drude", "drudeNbtholeCut", "Nonbonded Thole interactions "
       "interaction radius", &drudeNbtholeCut, 5.0);
   opts.range("drudeNbtholeCut", POSITIVE);
   opts.optionalB("drude", "drudeHardWall", "Apply maximum Drude bond length "
       "restriction?", &drudeHardWallOn, TRUE);
   opts.optional("drude", "drudeBondLen", "Drude oscillator bond length "
       "beyond which to apply restraint", &drudeBondLen, 0.25);
   opts.range("drudeBondLen", POSITIVE);
   opts.optional("drude", "drudeBondConst", "Drude oscillator restraining "
       "force constant", &drudeBondConst, 40000.0);
   opts.range("drudeBondConst", POSITIVE);

   // Pair interaction calculations
    opts.optionalB("main", "pairInteraction",
	"Are pair interactions calculated?", &pairInteractionOn, FALSE);
    opts.optional("pairInteraction", "pairInteractionFile",
	"PDB files with interaction flags " "default is the input PDB file",
	PARSE_STRING);
    opts.optional("pairInteraction", "pairInteractionCol",
	"Column in the pairInteractionFile with the interaction flags",
	PARSE_STRING);
    opts.require("pairInteraction", "pairInteractionGroup1",
        "Flag for interaction group 1", &pairInteractionGroup1);
    opts.optional("pairInteraction", "pairInteractionGroup2",
        "Flag for interaction group 2", &pairInteractionGroup2, -1);
    opts.optionalB("pairInteraction", "pairInteractionSelf",
        "Compute only within-group interactions?", &pairInteractionSelf,
        FALSE);
   // Options for CG simulations
   opts.optionalB("main", "cosAngles", "Are some angles cosine-based?", &cosAngles, FALSE);


   //  Dihedral angle dynamics
   opts.optionalB("main", "globalTest", "Should global integration (for development) be used?",
    &globalOn, FALSE);
   opts.optionalB("main", "dihedral", "Should dihedral angle dynamics be performed?",
    &dihedralOn, FALSE);
   COLDOn = FALSE;
   opts.optionalB("dihedral", "COLD", "Should overdamped Langevin dynamics be performed?",
    &COLDOn, FALSE);
   opts.require("COLD", "COLDTemp", "Temperature for heat bath in COLD",
    &COLDTemp);
   opts.range("COLDTemp", NOT_NEGATIVE);
   opts.units("COLDTemp", N_KELVIN);
   opts.require("COLD", "COLDRate", "Damping rate for COLD",
    &COLDRate, 3000.0);
   opts.range("COLDRate", NOT_NEGATIVE);

   //  Get the parameters for temperature coupling
   opts.optionalB("main", "tcouple",
      "Should temperature coupling be performed?",
      &tCoupleOn, FALSE);
   opts.require("tcouple", "tCoupleTemp",
    "Temperature for temperature coupling", &tCoupleTemp);
   opts.range("tCoupleTemp", NOT_NEGATIVE);
   opts.units("tCoupleTemp", N_KELVIN);
   opts.optional("tCouple", "tCoupleFile", "PDB file with temperature "
     "coupling terms (B(i)) (default is the PDB input file)",
     PARSE_STRING);
   opts.optional("tCouple", "tCoupleCol", "Column in the tCoupleFile "
     "containing the temperature coupling term B(i);\n"
     "default is 'O'", PARSE_STRING);

   opts.optionalB("main", "stochRescale",
      "Should stochastic velocity rescaling be performed?",
      &stochRescaleOn, FALSE);
   opts.require("stochRescale", "stochRescaleTemp",
      "Temperature for stochastic velocity rescaling",
       &stochRescaleTemp);
   opts.range("stochRescaleTemp", NOT_NEGATIVE);
   opts.units("stochRescaleTemp", N_KELVIN);
   opts.require("stochRescale", "stochRescalePeriod",
      "Time scale for stochastic velocity rescaling (ps)",
       &stochRescalePeriod);
   opts.range("stochRescalePeriod", POSITIVE);
   opts.optional("stochRescale", "stochRescaleFreq",
       "Number of steps between stochastic rescalings",
        &stochRescaleFreq);
   opts.range("stochRescaleFreq", POSITIVE);
   opts.optionalB("stochRescale", "stochRescaleHeat",
       "Should heat transfer and work be computed?", &stochRescaleHeat, FALSE);

   opts.optional("main", "rescaleFreq", "Number of steps between "
    "velocity rescaling", &rescaleFreq);
   opts.range("rescaleFreq", POSITIVE);
   opts.optional("main", "rescaleTemp", "Target temperature for velocity rescaling",
    &rescaleTemp);
   opts.range("rescaleTemp", NOT_NEGATIVE);
   opts.units("rescaleTemp", N_KELVIN);

   opts.optional("main", "reassignFreq", "Number of steps between "
    "velocity reassignment", &reassignFreq);
   opts.range("reassignFreq", POSITIVE);
   opts.optional("main", "reassignTemp", "Target temperature for velocity reassignment",
    &reassignTemp);
   opts.range("reassignTemp", NOT_NEGATIVE);
   opts.units("reassignTemp", N_KELVIN);
   opts.optional("main", "reassignIncr", "Temperature increment for velocity reassignment",
    &reassignIncr);
   opts.units("reassignIncr", N_KELVIN);
   opts.optional("main", "reassignHold", "Final holding temperature for velocity reassignment",
    &reassignHold);
   opts.range("reassignHold", NOT_NEGATIVE);
   opts.units("reassignHold", N_KELVIN);

   ////  Group rather than atomic pressure
   opts.optionalB("main", "useGroupPressure",
      "Use group rather than atomic quantities for pressure control?",
      &useGroupPressure, FALSE);

   ////  Anisotropic cell fluctuations
   opts.optionalB("main", "useFlexibleCell",
      "Use anisotropic cell fluctuation for pressure control?",
      &useFlexibleCell, FALSE);

   // Fix specific cell dimensions
   opts.optionalB("main", "fixCellDims",
      "Fix some cell dimensions?",
      &fixCellDims, FALSE);

   opts.optionalB("fixCellDims", "fixCellDimX",
      "Fix the X dimension?",
      &fixCellDimX, FALSE);
   opts.optionalB("fixCellDims", "fixCellDimY",
      "Fix the Y dimension?",
      &fixCellDimY, FALSE);
   opts.optionalB("fixCellDims", "fixCellDimZ",
      "Fix the Z dimension?",
      &fixCellDimZ, FALSE);

   ////  Constant dimension ratio in X-Y plane
   opts.optionalB("main", "useConstantRatio",
      "Use constant X-Y ratio for pressure control?",
      &useConstantRatio, FALSE);

   ////  Constant area and normal pressure conditions
   opts.optionalB("main", "useConstantArea",
      "Use constant area for pressure control?",
      &useConstantArea, FALSE);

   //// Exclude atoms from pressure
   opts.optionalB("main", "excludeFromPressure",
	"Should some atoms be excluded from pressure rescaling?",
	&excludeFromPressure, FALSE);
   opts.optional("excludeFromPressure", "excludeFromPressureFile",
	"PDB file for atoms to be excluded from pressure",
        PARSE_STRING);
   opts.optional("excludeFromPressure", "excludeFromPressureCol",
        "Column in the excludeFromPressureFile"
        "containing the flags (nonzero means excluded);\n"
        "default is 'O'", PARSE_STRING);

   ////  Berendsen pressure bath coupling
   opts.optionalB("main", "BerendsenPressure",
      "Should Berendsen pressure bath coupling be performed?",
      &berendsenPressureOn, FALSE);
   opts.require("BerendsenPressure", "BerendsenPressureTarget",
    "Target pressure for pressure coupling",
    &berendsenPressureTarget);
   // opts.units("BerendsenPressureTarget",);
   opts.require("BerendsenPressure", "BerendsenPressureCompressibility",
    "Isothermal compressibility for pressure coupling",
    &berendsenPressureCompressibility);
   // opts.units("BerendsenPressureCompressibility",);
   opts.require("BerendsenPressure", "BerendsenPressureRelaxationTime",
    "Relaxation time for pressure coupling",
    &berendsenPressureRelaxationTime);
   opts.range("BerendsenPressureRelaxationTime", POSITIVE);
   opts.units("BerendsenPressureRelaxationTime", N_FSEC);
   opts.optional("BerendsenPressure", "BerendsenPressureFreq",
    "Number of steps between volume rescaling",
    &berendsenPressureFreq, 1);
   opts.range("BerendsenPressureFreq", POSITIVE);

   ////  Langevin Piston pressure control
   opts.optionalB("main", "LangevinPiston",
      "Should Langevin piston pressure control be used?",
      &langevinPistonOn, FALSE);
   opts.optionalB("LangevinPiston", "LangevinPistonBarrier",
      "Should Langevin piston barrier be used?",
      &langevinPistonBarrier, TRUE);
   opts.require("LangevinPiston", "LangevinPistonTarget",
      "Target pressure for pressure control",
      &langevinPistonTarget);
   opts.require("LangevinPiston", "LangevinPistonPeriod",
      "Oscillation period for pressure control",
      &langevinPistonPeriod);
   opts.range("LangevinPistonPeriod", POSITIVE);
   opts.units("LangevinPistonPeriod", N_FSEC);
   opts.require("LangevinPiston", "LangevinPistonDecay",
      "Decay time for pressure control",
      &langevinPistonDecay);
   opts.range("LangevinPistonDecay", POSITIVE);
   opts.units("LangevinPistonDecay", N_FSEC);
   opts.require("LangevinPiston", "LangevinPistonTemp",
      "Temperature for pressure control piston",
      &langevinPistonTemp);
   opts.range("LangevinPistonTemp", POSITIVE);
   opts.units("LangevinPistonTemp", N_KELVIN);
   opts.optional("LangevinPiston", "StrainRate",
      "Initial strain rate for pressure control (x y z)",
      &strainRate);

   // Multigrator temperature and/or pressure control
   opts.optionalB("main", "Multigrator",
      "Should multigrator temperature and/or pressure control be used?",
      &multigratorOn, FALSE);
   opts.require("Multigrator", "MultigratorPressureTarget",
    "Target pressure for pressure coupling",
    &multigratorPressureTarget);
   opts.require("Multigrator", "MultigratorTemperatureTarget",
    "Target temperature for temperature coupling",
    &multigratorTemperatureTarget);
   opts.require("Multigrator", "MultigratorPressureFreq",
    "Number of steps between pressure control moves",
    &multigratorPressureFreq);
   opts.range("MultigratorPressureFreq", POSITIVE);
   opts.optional("Multigrator", "MultigratorPressureRelaxationTime",
    "Relaxation time for pressure coupling is fs",
    &multigratorPressureRelaxationTime, 30000);
   opts.range("MultigratorPressureRelaxationTime", POSITIVE);
   opts.units("MultigratorPressureRelaxationTime", N_FSEC);
   opts.optional("Multigrator", "MultigratorTemperatureRelaxationTime",
    "Relaxation time for temperature coupling is fs",
    &multigratorTemperatureRelaxationTime, 1000);
   opts.range("MultigratorTemperatureRelaxationTime", POSITIVE);
   opts.units("MultigratorTemperatureRelaxationTime", N_FSEC);
   opts.require("Multigrator", "MultigratorTemperatureFreq",
    "Number of steps between temperature control moves",
    &multigratorTemperatureFreq);
   opts.range("MultigratorTemperatureFreq", POSITIVE);
   opts.optional("Multigrator", "MultigratorNoseHooverChainLength",
    "Nose-Hoover chain length",
    &multigratorNoseHooverChainLength, 4);
   opts.range("MultigratorNoseHooverChainLength", POSITIVE);

   //// Surface tension
   opts.optional("main", "SurfaceTensionTarget",
      "Surface tension in the x-y plane",
      &surfaceTensionTarget, 0);

   //// Pressure Profile calculations
   opts.optionalB("main", "pressureprofile", "Compute pressure profile?",
     &pressureProfileOn, FALSE);
   opts.require("pressureprofile", "pressureprofileslabs",
     "Number of pressure profile slabs", &pressureProfileSlabs, 10);
   opts.optional("pressureprofile", "pressureprofilefreq",
     "How often to store profile data", &pressureProfileFreq, 1);
   opts.optional("pressureprofile", "pressureProfileAtomTypes",
     "Number of pressure profile atom types", &pressureProfileAtomTypes, 1);
   opts.range("pressureProfileAtomTypes", POSITIVE);
   opts.optional("pressureProfile", "pressureProfileAtomTypesFile",
	"PDB files with pressure profile atom types" "default is the input PDB file",
	PARSE_STRING);
   opts.optional("pressureProfile", "pressureProfileAtomTypesCol",
	"Column in the pressureProfileAtomTypesFile with the atom types ",
	PARSE_STRING);
   opts.optionalB("pressureProfile", "pressureProfileEwald",
       "Compute Ewald contribution to pressure profile",
       &pressureProfileEwaldOn, FALSE);
   opts.optional("pressureProfile", "pressureProfileEwaldX",
       "Ewald grid size X", &pressureProfileEwaldX, 10);
   opts.range("pressureProfileEwaldX", POSITIVE);
   opts.optional("pressureProfile", "pressureProfileEwaldY",
       "Ewald grid size Y", &pressureProfileEwaldY, 10);
   opts.range("pressureProfileEwaldY", POSITIVE);
   opts.optional("pressureProfile", "pressureProfileEwaldZ",
       "Ewald grid size Z", &pressureProfileEwaldZ, 10);
   opts.range("pressureProfileEwaldZ", POSITIVE);

   /// accelerated MD parameters
   opts.optionalB("main", "accelMD", "Perform acclerated MD?", &accelMDOn, FALSE);
   opts.optional("accelMD", "accelMDFirstStep", "First accelMD step", &accelMDFirstStep, 0);
   opts.range("accelMDFirstStep", NOT_NEGATIVE);
   opts.optional("accelMD", "accelMDLastStep", "Last accelMD step", &accelMDLastStep, 0);
   opts.range("accelMDLastStep", NOT_NEGATIVE);
   opts.optional("accelMD", "accelMDOutFreq", "Frequency of accelMD output", &accelMDOutFreq, 1);
   opts.range("accelMDOutFreq", POSITIVE);
   opts.optionalB("accelMD", "accelMDdihe", "Apply boost to dihedral potential", &accelMDdihe, TRUE);
   opts.optionalB("accelMD", "accelMDDebugOn", "Debugging accelMD", &accelMDDebugOn, FALSE);
   opts.optional("accelMD", "accelMDE","E for AMD", &accelMDE);
   opts.units("accelMDE", N_KCAL);
   opts.optional("accelMD", "accelMDalpha","alpha for AMD", &accelMDalpha);
   opts.units("accelMDalpha", N_KCAL);
   opts.range("accelMDalpha", POSITIVE);
   opts.optionalB("accelMD", "accelMDdual", "Apply dual boost", &accelMDdual, FALSE);
   opts.optional("accelMDdual", "accelMDTE","E for total potential under accelMDdual mode", &accelMDTE);
   opts.units("accelMDTE", N_KCAL);
   opts.optional("accelMDdual", "accelMDTalpha","alpha for total potential under accelMDdual mode", &accelMDTalpha);
   opts.units("accelMDTalpha", N_KCAL);
   opts.range("accelMDTalpha", POSITIVE);
   // GaMD parameters
   opts.optionalB("accelMD", "accelMDG", "Perform Gaussian accelMD calculation?", &accelMDG, FALSE);
   opts.optional("accelMDG", "accelMDGiE", "Flag to set the mode iE in Gaussian accelMD", &accelMDGiE, 1);
   opts.optional("accelMDG", "accelMDGcMDSteps", "Number of cMD steps", &accelMDGcMDSteps, 1000000);
   opts.range("accelMDGcMDSteps", NOT_NEGATIVE);
   opts.optional("accelMDG", "accelMDGEquiSteps", "Number of equilibration steps after adding boost potential", &accelMDGEquiSteps, 1000000);
   opts.range("accelMDGEquiSteps", NOT_NEGATIVE);
   opts.require("accelMDG", "accelMDGcMDPrepSteps", "Number of preparation cMD steps", &accelMDGcMDPrepSteps, 200000);
   opts.range("accelMDGcMDPrepSteps", NOT_NEGATIVE);
   opts.require("accelMDG", "accelMDGEquiPrepSteps", "Number of preparation equilibration steps", &accelMDGEquiPrepSteps, 200000);
   opts.range("accelMDGEquiPrepSteps", NOT_NEGATIVE);
   opts.optional("accelMDG", "accelMDGStatWindow", "Number of steps to calculate avg and std", &accelMDGStatWindow, -1);
   opts.optional("accelMDG", "accelMDGSigma0P", "Upper limit of std of total potential", &accelMDGSigma0P, 6.0);
   opts.units("accelMDGSigma0P", N_KCAL);
   opts.range("accelMDGSigma0P", NOT_NEGATIVE);
   opts.optional("accelMDG", "accelMDGSigma0D", "Upper limit of std of dihedral potential", &accelMDGSigma0D, 6.0);
   opts.units("accelMDGSigma0D", N_KCAL);
   opts.range("accelMDGSigma0D", NOT_NEGATIVE);
   opts.optionalB("accelMDG", "accelMDGRestart", "Flag to set use restart file in Gaussian accelMD", &accelMDGRestart, FALSE);
   opts.require("accelMDGRestart", "accelMDGRestartFile", "Restart file name for Gaussian accelMD", accelMDGRestartFile);
   opts.optionalB("accelMDG", "accelMDGresetVaftercmd", "Flag to reset potential after accelMDGcMDSteps steps",
	   &accelMDGresetVaftercmd, FALSE);

   // Adaptive Temperature Sampling (adaptTemp) parameters
   opts.optionalB("main", "adaptTempMD", "Perform adaptive temperature sampling", &adaptTempOn, FALSE);
   opts.optional("adaptTempMD", "adaptTempFirstStep", "First adaptTemp step", &adaptTempFirstStep, 0);
   opts.range("adaptTempFirstStep", NOT_NEGATIVE);
   opts.optional("adaptTempMD", "adaptTempLastStep", "Last adaptTemp step", &adaptTempLastStep, 0);
   opts.range("adaptTempLastStep", NOT_NEGATIVE);
   opts.optional("adaptTempMD", "adaptTempOutFreq", "Frequency of adaptTemp output", &adaptTempOutFreq, 10);
   opts.range("adaptTempOutFreq", POSITIVE);
   opts.optional("adaptTempMD", "adaptTempFreq", "Frequency of writing average energies to adaptTempOutFile", &adaptTempFreq, 10);
   opts.range("adaptTempFreq", POSITIVE);
   opts.optionalB("adaptTempMD", "adaptTempDebug", "Print debug output for adaptTemp", &adaptTempDebug, FALSE);
   opts.optional("adaptTempMD", "adaptTempTmin","Minimun temperature for adaptTemp", &adaptTempTmin);
   opts.units("adaptTempTmin", N_KELVIN);
   opts.range("adaptTempTmin", POSITIVE);
   opts.optional("adaptTempMD", "adaptTempTmax","Maximum temperature for adaptTemp", &adaptTempTmax);
   opts.units("adaptTempTmax", N_KELVIN);
   opts.range("adaptTempTmax", POSITIVE);
   opts.optional("adaptTempMD", "adaptTempBins","Number of bins to store average energies", &adaptTempBins,0);
   opts.range("adaptTempBins", NOT_NEGATIVE);
   opts.optional("adaptTempMD", "adaptTempDt", "Integration timestep for Temp. updates", &adaptTempDt, 0.0001);
   opts.units("adaptTempDt", N_FSEC);
   opts.range("adaptTempDt", NOT_NEGATIVE);
   opts.optional("adaptTempMD", "adaptTempAutoDt", "Average temperature update in percent of temperature range", &adaptTempAutoDt, 0.0);
   opts.range("adaptTempAutoDt", NOT_NEGATIVE);
   opts.optional("adaptTempMD", "adaptTempCgamma", "Adaptive bin averaging constant", &adaptTempCgamma, 0.1);
   opts.range("adaptTempCgamma", NOT_NEGATIVE);
   opts.optionalB("adaptTempMD","adaptTempLangevin","Send adaptTemp temperature to langevin thermostat",&adaptTempLangevin,TRUE);
   opts.optionalB("adaptTempMD","adaptTempRescaling","Send adaptTemp temperature to velocity rescaling thermostat", &adaptTempRescale,TRUE);
   opts.optional("adaptTempMD", "adaptTempInFile", "File containing restart information for adaptTemp", adaptTempInFile);
   opts.optional("adaptTempMD", "adaptTempRestartFile", "File for writing adaptTemp restart information", adaptTempRestartFile);
   opts.require("adaptTempRestartFile","adaptTempRestartFreq", "Frequency of writing restart file", &adaptTempRestartFreq,0);
   opts.range("adaptTempRestartFreq",NOT_NEGATIVE);
   opts.optionalB("adaptTempMD", "adaptTempRandom", "Randomly assign a temperature if we step out of range", &adaptTempRandom, FALSE);
}

void SimParameters::config_parser_constraints(ParseOptions &opts) {

   ////  Fixed Atoms
   opts.optionalB("main", "fixedatoms", "Are there fixed atoms?",
    &fixedAtomsOn, FALSE);
   opts.optionalB("fixedatoms", "fixedAtomsForces",
     "Calculate forces between fixed atoms?  (Required to unfix during run.)",
     &fixedAtomsForces, FALSE);
   opts.optional("fixedatoms", "fixedAtomsFile", "PDB file with flags for "
     "fixed atoms (default is the PDB input file)",
     PARSE_STRING);
   opts.optional("fixedatoms", "fixedAtomsCol", "Column in the fixedAtomsFile "
     "containing the flags (nonzero means fixed);\n"
     "default is 'O'", PARSE_STRING);
   opts.optional("fixedatoms", "fixedAtomListFile", "the text input file for fixed atoms "
                 "used for parallel input IO", PARSE_STRING);
   opts.optionalB("fixedatoms", "fixedAtomsForceOutput",
     "Do we write out forces acting on fixed atoms?",
     &fixedAtomsForceOutput, FALSE);

   ////  Harmonic Constraints
   opts.optionalB("main", "constraints", "Are harmonic constraints active?",
     &constraintsOn, FALSE);
   opts.require("constraints", "consexp", "Exponent for harmonic potential",
    &constraintExp, 2);
   opts.range("consexp", POSITIVE);
#ifndef MEM_OPT_VERSION
   opts.require("constraints", "consref", "PDB file containing reference "
    "positions",
    PARSE_STRING);
   opts.require("constraints", "conskfile", "PDB file containing force "
    "constaints in one of the columns", PARSE_STRING);
   opts.require("constraints", "conskcol", "Column of conskfile to use "
    "for the force constants", PARSE_STRING);
#else
   opts.require("constraints", "consAtomListFile", "the text input file for constrained atoms "
                 "used for parallel input IO", PARSE_STRING);
#endif
   opts.require("constraints", "constraintScaling", "constraint scaling factor",
     &constraintScaling, 1.0);
   opts.range("constraintScaling", NOT_NEGATIVE);



   //****** BEGIN selective restraints (X,Y,Z) changes

   //// selective restraints (X,Y,Z)
   opts.optionalB("constraints", "selectConstraints",
   "Restrain only selected Cartesian components of the coordinates?",
     &selectConstraintsOn, FALSE);
   opts.optionalB("selectConstraints", "selectConstrX",
   "Restrain X components of coordinates ", &constrXOn, FALSE);
   opts.optionalB("selectConstraints", "selectConstrY",
   "Restrain Y components of coordinates ", &constrYOn, FALSE);
   opts.optionalB("selectConstraints", "selectConstrZ",
   "Restrain Z components of coordinates ", &constrZOn, FALSE);
   //****** END selective restraints (X,Y,Z) changes

   // spherical constraints
   opts.optionalB("constraints", "sphericalConstraints",
   "Restrain only radial spherical component of the coordinates?",
     &sphericalConstraintsOn, FALSE);
   opts.optional("sphericalConstraints", "sphericalConstrCenter",
   "Center of spherical constraints", &sphericalConstrCenter);

   //****** BEGIN moving constraints changes

   //// Moving Harmonic Constraints
   opts.optionalB("constraints", "movingConstraints",
      "Are some of the constraints moving?",
      &movingConstraintsOn, FALSE);
   opts.require("movingConstraints", "movingConsVel",
    "Velocity of the movement, A/timestep", &movingConsVel);
   //****** END moving constraints changes

   // BEGIN rotating constraints changes
   opts.optionalB("constraints", "rotConstraints",
      "Are the constraints rotating?",
      &rotConstraintsOn, FALSE);
   opts.require("rotConstraints", "rotConsAxis",
    "Axis of rotation", &rotConsAxis);
   opts.require("rotConstraints", "rotConsPivot",
    "Pivot point of rotation",
    &rotConsPivot);
   opts.require("rotConstraints", "rotConsVel",
    "Velocity of rotation, deg/timestep", &rotConsVel);

   // END rotating constraints changes

   // external command forces
   opts.optionalB("main", "extForces", "External command forces?",
      &extForcesOn, FALSE);
   opts.require("extForces", "extForcesCommand",
      "External forces command", extForcesCommand);
   opts.require("extForces", "extCoordFilename",
      "External forces coordinate filename", extCoordFilename);
   opts.require("extForces", "extForceFilename",
      "External forces force filename", extForceFilename);


  // QM/MM forces
   opts.optionalB("main", "QMForces", "Apply QM forces?",
      &qmForcesOn, FALSE);
   opts.require("QMForces", "QMSoftware",
      "software whose format will be used for input/output", qmSoftware);
   opts.require("QMForces", "QMExecPath",
      "path to executable", qmExecPath);
   opts.optional("QMForces", "QMChargeMode",
      "type of QM atom charges gathered from the QM software", qmChrgModeS);
   opts.require("QMForces", "QMColumn",
      "column defining QM and MM regions", qmColumn);
   opts.require("QMForces", "QMBaseDir",
      "base path and name for QM input and output (preferably in memory)", qmBaseDir);
   opts.optional("QMForces", "QMConfigLine",
      "Configuration line for QM (multiple inputs allowed)", PARSE_MULTIPLES);
   opts.optional("QMForces", "QMParamPDB",
      "PDB with QM parameters", qmParamPDB);
   opts.optional("QMForces", "QMPrepProc",
      "initial preparation executable", qmPrepProc);
   opts.optional("QMForces", "QMSecProc",
      "secondary executable", qmSecProc);
   opts.optional("QMForces", "QMCharge",
      "charge of the QM group", PARSE_MULTIPLES);
   opts.optionalB("QMForces", "QMChargeFromPSF",
      "gets charge of the QM group form PSF values", &qmChrgFromPSF, FALSE);
   opts.optional("QMForces", "QMMult",
      "multiplicity of the QM group", PARSE_MULTIPLES);
   opts.optional("QMForces", "QMLinkElement",
      "element of link atom", PARSE_MULTIPLES);
   opts.optionalB("QMForces", "QMReplaceAll",
      "replace all NAMD forces with QM forces", &qmReplaceAll, FALSE);
   opts.optional("QMForces", "QMPCStride",
      "frequency of selection of point charges", &qmPCSelFreq, 1);
   opts.range("QMPCStride", POSITIVE);
   opts.optionalB("QMForces", "QMNoPntChrg",
      "no point charges will be passed to the QM system(s)", &qmNoPC, FALSE);
   opts.optionalB("QMForces", "QMElecEmbed",
      "activates electrostatic embedding", &qmElecEmbed, TRUE);
   opts.optionalB("QMForces", "QMVdWParams",
      "use special VdW parameters for QM atoms", &qmVDW, FALSE);
   opts.optional("QMForces", "QMBondColumn",
      "column defining QM-MM bomnds", qmBondColumn);
   opts.optionalB("QMForces", "QMBondDist",
      "values in QMBondColumn defines the distance of new link atom", &qmBondDist, FALSE);
   opts.optional("QMForces", "QMBondValueType",
      "type of value in bond column: len or ratio", qmBondValueTypeS);
   opts.optional("QMForces", "QMBondScheme",
      "type of treatment given to QM-MM bonds.", qmBondSchemeS);
   opts.optional("QMForces", "QMenergyStride",
      "frequency of QM specific energy output (every x steps)", &qmEnergyOutFreq, 1);
   opts.optional("QMForces", "QMOutStride",
      "frequency of QM specific charge output (every x steps)", &qmOutFreq, 0);
   opts.range("QMOutStride", NOT_NEGATIVE);
   opts.optional("QMForces", "QMPositionOutStride",
      "frequency of QM specific position output (every x steps)", &qmPosOutFreq, 0);
   opts.range("QMPositionOutStride", NOT_NEGATIVE);
   opts.optional("QMForces", "QMSimsPerNode",
      "QM executions per node", &qmSimsPerNode, 1);
   opts.range("QMSimsPerNode", POSITIVE);
   opts.optionalB("QMForces", "QMSwitching",
      "apply switching to point charges.", &qmPCSwitchOn, FALSE);
   opts.optional("QMForces", "QMSwitchingType",
      "How are charges scaled down to be presented to QM groups.", qmPCSwitchTypeS);
   opts.optional("QMForces", "QMPointChargeScheme",
      "type of treatment given to the total sum of point charges.", qmPCSchemeS);
   opts.optionalB("QMForces", "QMCustomPCSelection",
      "custom and fixed selection of point charges per QM group.", &qmCustomPCSel, FALSE);
   opts.optional("QMForces", "QMCustomPCFile",
      "file with a selection of point charges for a single QM group", PARSE_MULTIPLES);
   opts.optionalB("QMForces", "QMLiveSolventSel",
      "Continuously update the selection of solvent molecules in QM groups", &qmLSSOn, FALSE);
   opts.optional("QMForces", "QMLSSFreq",
      "frequency of QM water selection update", &qmLSSFreq, 100);
   opts.range("QMLSSFreq", POSITIVE);
   opts.optional("QMForces", "QMLSSResname",
      "residue name for the solvent molecules (TIP3).", qmLSSResname);
   opts.optional("QMForces", "QMLSSMode",
      "mode of selection of point solvent molecules", qmLSSModeS);
   opts.optional("QMForces", "QMLSSRef",
      "for COM mode, defines reference for COM distance calculation", PARSE_MULTIPLES);
   opts.optionalB("QMForces", "QMCSMD",
      "Do we use Conditional SMD option?", &qmCSMD, FALSE);
   opts.optional("QMForces", "QMCSMDFile",
                "File for Conditional SMD information",qmCSMDFile);

   //print which bad contacts are being moved downhill
   opts.optionalB("main", "printBadContacts", "Print atoms with huge forces?",
      &printBadContacts, FALSE);

   /* GBIS generalized born implicit solvent*/

   opts.optionalB("main", "GBIS", "Use GB implicit solvent?",
      &GBISOn, FALSE);
   opts.optionalB("main", "GBISSer", "Use GB implicit solvent?",
      &GBISserOn, FALSE);

   opts.optional("GBIS", "solventDielectric",
      "Solvent Dielectric", &solvent_dielectric, 78.5);
   opts.optional("GBIS", "intrinsicRadiusOffset",
      "Coulomb Radius Offset", &coulomb_radius_offset, 0.09);
   opts.optional("GBIS", "ionConcentration",
      "Ion Concentration", &ion_concentration, 0.2); //0.2 mol/L
   opts.optional("GBIS", "GBISDelta",
      "delta from GBOBC", &gbis_delta, 1.0); //0.8 or 1.0
   opts.optional("GBIS", "GBISBeta",
      "beta from GBOBC", &gbis_beta, 0.8);   //0.0 or 0.8
   opts.optional("GBIS", "GBISGamma",
      "gamma from GBOBC", &gbis_gamma, 4.85);//2.290912 or 4.85
   opts.optional("GBIS", "alphaCutoff",
      "cutoff for calculating effective born radius", &alpha_cutoff, 15);
   opts.optional("GBIS", "alphaMax",
      "maximum allowable born radius", &alpha_max, 30);
   opts.optional("GBIS", "fsMax",
      "maximum screened intrinsic radius", &fsMax, 1.728);

   opts.optionalB("main", "SASA", "Use Linear Combination of Pairwise Overlaps (LCPO) for calculating SASA",
      &LCPOOn, FALSE);
   opts.optional("SASA", "surfaceTension",
      "Surfce Tension for SASA (kcal/mol/Ang^2)", &surface_tension, 0.005);

   //****** BEGIN SMD constraints changes

   // SMD constraints
   opts.optionalB("main", "SMD",
      "Do we use SMD option?",
      &SMDOn, FALSE);
   opts.require("SMD", "SMDVel",
		"Velocity of the movement, A/timestep", &SMDVel);
   opts.range("SMDVel", NOT_NEGATIVE);
   opts.require("SMD", "SMDDir",
		"Direction of movement", &SMDDir);
   opts.require("SMD", "SMDk",
                "Elastic constant for SMD", &SMDk);
   opts.optional("SMD", "SMDk2",
                "Transverse elastic constant for SMD", &SMDk2, 0);
   opts.range("SMDk", NOT_NEGATIVE);
   opts.range("SMDk2", NOT_NEGATIVE);
   opts.require("SMD", "SMDFile",
		"File for SMD information",
                 SMDFile);
   opts.optional("SMD", "SMDOutputFreq",
		 "Frequency of output",
		 &SMDOutputFreq, 1);
   opts.range("SMDOutputFreq", POSITIVE);

   //****** END SMD constraints changes

   //****** BEGIN tabulated energies section
   opts.optionalB("main", "tabulatedEnergies", "Do we get energies from a table?", &tabulatedEnergies, FALSE);
//   opts.require("tabulatedEnergies", "tableNumTypes","Number of types for energy tabulation", &tableNumTypes);
   opts.require("tabulatedEnergies", "tabulatedEnergiesFile", "File containing energy table", tabulatedEnergiesFile);
   opts.require("tabulatedEnergies", "tableInterpType", "Cubic or linear interpolation", tableInterpType);

   // TMD parameters
   opts.optionalB("main", "TMD", "Perform Targeted MD?", &TMDOn, FALSE);
   opts.optional("TMD", "TMDk", "Elastic constant for TMD", &TMDk, 0);
   opts.range("TMDk", NOT_NEGATIVE);
   opts.require("TMD", "TMDFile", "File for TMD information", TMDFile);
   opts.optionalB("TMD", "TMDDiffRMSD", "Restrain Difference between the RMSD from two structures", &TMDDiffRMSD, FALSE);
   opts.require("TMDDiffRMSD", "TMDFile2",  "Second file for TMD information", TMDFile2);

   opts.optional("TMD", "TMDOutputFreq", "Frequency of TMD output",
       &TMDOutputFreq, 1);
   opts.range("TMDOutputFreq", POSITIVE);
   opts.require("TMD", "TMDLastStep", "Last TMD timestep", &TMDLastStep);
   opts.range("TMDLastStep", POSITIVE);
   opts.optional("TMD", "TMDFirstStep", "First TMD step (default 0)", &TMDFirstStep, 0);
   opts.optional("TMD", "TMDInitialRMSD", "Target RMSD at first TMD step (default -1 to use initial coordinates)", &TMDInitialRMSD);
   TMDInitialRMSD = -1;
   opts.optional("TMD", "TMDFinalRMSD", "Target RMSD at last TMD step (default 0 )", &TMDFinalRMSD, 0);
   opts.range("TMDInitialRMSD", NOT_NEGATIVE);
   // End of TMD parameters

   // Symmetry restraint parameters
   opts.optionalB("main", "symmetryRestraints", "Enable symmetry restraints?", &symmetryOn, FALSE);
   opts.optional("symmetryRestraints", "symmetryk", "Elastic constant for symmetry restraints", &symmetryk, 0);
   opts.range("symmetryk", NOT_NEGATIVE);
   opts.optional("symmetryRestraints", "symmetrykfile", "PDB file specifying force contants on a per-atom basis", PARSE_MULTIPLES);
   opts.optionalB("symmetryRestraints", "symmetryScaleForces", "Scale applied forces over time?", &symmetryScaleForces, FALSE);
   opts.require("symmetryRestraints", "symmetryFile", "File for symmetry information", PARSE_MULTIPLES);
   opts.optional("symmetryRestraints", "symmetryMatrixFile", "File(s) for transfromation matrices", PARSE_MULTIPLES);
   opts.optional("symmetryRestraints", "symmetryLastStep", "Last symmetry timestep", &symmetryLastStep, -1);
   opts.optional("symmetryRestraints", "symmetryFirstStep", "First symmetry step (default 0)", &symmetryFirstStep, 0);
   opts.optional("symmetryRestraints", "symmetryLastFullStep", "Last full force symmetry timestep (default symmetryLastStep)", &symmetryLastFullStep, symmetryLastStep);
   opts.optional("symmetryRestraints", "symmetryFirstFullStep", "First full force symmetry step (default symmetryFirstStep)", &symmetryFirstFullStep, symmetryFirstStep);
  //End of symmetry restraint parameters.

   ////  Global Forces / Tcl
   opts.optionalB("main", "tclForces", "Are Tcl global forces active?",
     &tclForcesOn, FALSE);
   opts.require("tclForces", "tclForcesScript",
     "Tcl script for global forces", PARSE_MULTIPLES);

   ////  Boundary Forces / Tcl
   opts.optionalB("main", "tclBC", "Are Tcl boundary forces active?",
     &tclBCOn, FALSE);
   opts.require("tclBC", "tclBCScript",
     "Tcl script defining calcforces for boundary forces", PARSE_STRING);
   tclBCScript = 0;
   opts.optional("tclBC", "tclBCArgs", "Extra args for calcforces command",
     tclBCArgs);
   tclBCArgs[0] = 0;

   ////  Global Forces / Misc
   opts.optionalB("main", "miscForces", "Are misc global forces active?",
     &miscForcesOn, FALSE);
   opts.optional("miscForces", "miscForcesScript",
     "script for misc forces", PARSE_MULTIPLES);

   ////  Free Energy Perturbation
   opts.optionalB("main", "freeEnergy", "Perform free energy perturbation?",
     &freeEnergyOn, FALSE);
   opts.require("freeEnergy", "freeEnergyConfig",
     "Configuration file for free energy perturbation", PARSE_MULTIPLES);

   ////  Constant Force
   opts.optionalB("main", "constantforce", "Apply constant force?",
     &consForceOn, FALSE);
   opts.optional("constantforce", "consForceFile",
       "Configuration file for constant forces", PARSE_STRING);
   opts.require("constantforce", "consForceScaling",
       "Scaling factor for constant forces", &consForceScaling, 1.0);

    //// Collective variables
    opts.optionalB("main", "colvars", "Is the colvars module enabled?",
      &colvarsOn, FALSE);
    opts.optional("colvars", "colvarsConfig",
      "configuration for the collective variables", PARSE_STRING);
    opts.optional("colvars", "colvarsInput",
      "input restart file for the collective variables", PARSE_STRING);

}

#ifdef OPENATOM_VERSION
void SimParameters::config_parser_openatom(ParseOptions &opts) {
  opts.optionalB("main", "openatom", "OpenAtom active?", &openatomOn, FALSE);
  opts.require("openatom", "openatomDriverFile", "What config file specifies openatom input parameters", PARSE_STRING);
  opts.require("openatom", "openatomPhysicsFile", "What structure file specifies openatom input system", PARSE_STRING);
  opts.require("openatom", "openatomPdbFile", "NAMD input file defining QM and MM regions", PARSE_STRING);
   opts.optional("openatom", "openatomCol", "Column in the openatomPdb with the QM/MM flag", PARSE_STRING);
}
#endif // OPENATOM_VERSION

/* BEGIN gf */
void SimParameters::config_parser_mgridforce(ParseOptions &opts) {
    //// Gridforce
    opts.optionalB("main", "mgridforce", "Is Multiple gridforce active?",
		   &mgridforceOn, FALSE);
    opts.optional("mgridforce", "mgridforcevolts", "Is Gridforce using Volts/eV as units?",
                  PARSE_MULTIPLES);
    opts.require("mgridforce", "mgridforcescale", "Scale factor by which to multiply "
		 "grid forces", PARSE_MULTIPLES);
    opts.require("mgridforce", "mgridforcefile", "PDB file containing force "
		 "multipliers in one of the columns", PARSE_MULTIPLES);
    opts.require("mgridforce", "mgridforcecol", "Column of gridforcefile to "
		 "use for force multiplier", PARSE_MULTIPLES);
    opts.optional("mgridforce", "mgridforcechargecol", "Column of gridforcefile to "
		  "use for charge", PARSE_MULTIPLES);
    opts.require("mgridforce", "mgridforcepotfile", "Gridforce potential file",
		 PARSE_MULTIPLES);
    opts.optional("mgridforce", "mgridforcecont1", "Use continuous grid "
		   "in K1 direction?", PARSE_MULTIPLES);
    opts.optional("mgridforce", "mgridforcecont2", "Use continuous grid "
		   "in K2 direction?", PARSE_MULTIPLES);
    opts.optional("mgridforce", "mgridforcecont3", "Use continuous grid "
		   "in K3 direction?", PARSE_MULTIPLES);
    opts.optional("mgridforce", "mgridforcevoff", "Gridforce potential offsets",
                  PARSE_MULTIPLES);
    opts.optional("mgridforce", "mgridforcelite", "Use Gridforce Lite?",
		  PARSE_MULTIPLES);
    opts.optional("mgridforce", "mgridforcechecksize", "Check if grid exceeds PBC cell dimensions?", PARSE_MULTIPLES);
}

void SimParameters::config_parser_gridforce(ParseOptions &opts) {
    //// Gridforce
    opts.optionalB("main", "gridforce", "Is Gridforce active?",
		   &gridforceOn, FALSE);
    opts.optionalB("gridforce", "gridforcevolts", "Is Gridforce using Volts/eV as units?",
		   &gridforceVolts, FALSE);
    opts.require("gridforce", "gridforcescale", "Scale factor by which to multiply "
		 "grid forces", &gridforceScale);
    opts.require("gridforce", "gridforcefile", "PDB file containing force "
		 "multipliers in one of the columns", PARSE_STRING);
    opts.require("gridforce", "gridforcecol", "Column of gridforcefile to "
		 "use for force multiplier", PARSE_STRING);
    opts.optional("gridforce", "gridforcechargecol", "Column of gridforcefile to "
		  "use for charge", PARSE_STRING);
    opts.require("gridforce", "gridforcepotfile", "Gridforce potential file",
		 PARSE_STRING);
    opts.optionalB("gridforce", "gridforcecont1", "Use continuous grid "
		   "in A1 direction?", &gridforceContA1, FALSE);
    opts.optionalB("gridforce", "gridforcecont2", "Use continuous grid "
		   "in A2 direction?", &gridforceContA2, FALSE);
    opts.optionalB("gridforce", "gridforcecont3", "Use continuous grid "
		   "in A3 direction?", &gridforceContA3, FALSE);
    opts.optional("gridforce", "gridforcevoff", "Gridforce potential offsets",
		  &gridforceVOffset);
    opts.optionalB("gridforce", "gridforcelite", "Use Gridforce Lite?",
		   &gridforceLite, FALSE);
    opts.optionalB("gridforce", "gridforcechecksize", "Check if grid exceeds PBC cell dimensions?",
		   &gridforcechecksize, TRUE);
}
/* END gf */

void SimParameters::config_parser_movdrag(ParseOptions &opts) {
   //// moving drag
   opts.optionalB("main", "movDragOn", "Do we apply moving drag?",
      &movDragOn, FALSE);
   opts.require("movDragOn", "movDragFile",
      "Main moving drag PDB file", movDragFile);
   opts.require("movDragOn", "movDragCol",
      "Main moving drag PDB column", PARSE_STRING);
   opts.require("movDragOn", "movDragGlobVel",
      "Global moving drag velocity (A/step)", &movDragGlobVel);
   opts.require("movDragOn", "movDragVelFile",
      "Moving drag linear velocity file", movDragVelFile);
}

void SimParameters::config_parser_rotdrag(ParseOptions &opts) {
   //// rotating drag
   opts.optionalB("main", "rotDragOn", "Do we apply rotating drag?",
      &rotDragOn, FALSE);
   opts.require("rotDragOn", "rotDragFile",
      "Main rotating drag PDB file", rotDragFile);
   opts.require("rotDragOn", "rotDragCol",
      "Main rotating drag PDB column", PARSE_STRING);
   opts.require("rotDragOn", "rotDragAxisFile",
      "Rotating drag axis file", rotDragAxisFile);
   opts.require("rotDragOn", "rotDragPivotFile",
      "Rotating drag pivot point file", rotDragPivotFile);
   opts.require("rotDragOn", "rotDragGlobVel",
      "Global rotating drag angular velocity (deg/step)", &rotDragGlobVel);
   opts.require("rotDragOn", "rotDragVelFile",
      "Rotating drag angular velocity file", rotDragVelFile);
   opts.require("rotDragOn", "rotDragVelCol",
      "Rotating drag angular velocity column", PARSE_STRING);
}

void SimParameters::config_parser_constorque(ParseOptions &opts) {
   //// "constant" torque
   opts.optionalB("main", "consTorqueOn", "Do we apply \"constant\" torque?",
      &consTorqueOn, FALSE);
   opts.require("consTorqueOn", "consTorqueFile",
      "Main \"constant\" torque PDB file", consTorqueFile);
   opts.require("consTorqueOn", "consTorqueCol",
      "Main \"constant\" torque PDB column", PARSE_STRING);
   opts.require("consTorqueOn", "consTorqueAxisFile",
      "\"Constant\" torque axis file", consTorqueAxisFile);
   opts.require("consTorqueOn", "consTorquePivotFile",
      "\"Constant\" torque pivot point file", consTorquePivotFile);
   opts.require("consTorqueOn", "consTorqueGlobVal",
      "Global \"constant\" torque value (Kcal/(mol*A^2))", &consTorqueGlobVal);
   opts.require("consTorqueOn", "consTorqueValFile",
      "\"constant\" torque factors file", consTorqueValFile);
   opts.require("consTorqueOn", "consTorqueValCol",
      "\"constant\" torque factors column", PARSE_STRING);
}

void SimParameters::config_parser_boundary(ParseOptions &opts) {

   //// Spherical Boundary Conditions
   opts.optionalB("main", "sphericalBC", "Are spherical boundary counditions "
      "active?", &sphericalBCOn, FALSE);
   opts.require("sphericalBC", "sphericalBCCenter",
     "Center of spherical boundaries", &sphericalCenter);
   opts.require("sphericalBC", "sphericalBCr1", "Radius for first sphere "
     "potential", &sphericalBCr1);
   opts.range("sphericalBCr1", POSITIVE);
   opts.units("sphericalBCr1", N_ANGSTROM);
   opts.require("sphericalBC", "sphericalBCk1", "Force constant for first "
    "sphere potential (+ is an inward force, - outward)",
    &sphericalBCk1);
   opts.units("sphericalBCk1", N_KCAL);
   opts.optional("sphericalBC", "sphericalBCexp1", "Exponent for first "
    "sphere potential", &sphericalBCexp1, 2);
   opts.range("sphericalBCexp1", POSITIVE);

   opts.optional("sphericalBCr1", "sphericalBCr2", "Radius for second sphere "
     "potential", &sphericalBCr2);
   opts.range("sphericalBCr2", POSITIVE);
   opts.units("sphericalBCr2", N_ANGSTROM);
   opts.require("sphericalBCr2", "sphericalBCk2", "Force constant for second "
    "sphere potential (+ is an inward force, - outward)",
    &sphericalBCk2);
   opts.units("sphericalBCk2", N_KCAL);
   opts.optional("sphericalBCr2", "sphericalBCexp2", "Exponent for second "
    "sphere potential", &sphericalBCexp2, 2);
   opts.range("sphericalBCexp2", POSITIVE);

   /////////////// Cylindrical Boundary Conditions
   opts.optionalB("main", "cylindricalBC", "Are cylindrical boundary counditions "
                  "active?", &cylindricalBCOn, FALSE);
   opts.require("cylindricalBC", "cylindricalBCr1", "Radius for first cylinder "
                 "potential", &cylindricalBCr1);
   opts.range("cylindricalBCr1", POSITIVE);
   opts.units("cylindricalBCr1", N_ANGSTROM);
   opts.require("cylindricalBC", "cylindricalBCk1", "Force constant for first "
                "cylinder potential (+ is an inward force, - outward)",
                &cylindricalBCk1);
   opts.units("cylindricalBCk1", N_KCAL);
   opts.optional("cylindricalBC", "cylindricalBCexp1", "Exponent for first "
                "cylinder potential", &cylindricalBCexp1, 2);
   opts.range("cylindricalBCexp1", POSITIVE);


// additions beyond those already found in spherical parameters    JJU
   opts.optional("cylindricalBC", "cylindricalBCAxis", "Cylinder axis (defaults to x)",
    PARSE_STRING);
   opts.require("cylindricalBC", "cylindricalBCCenter",
     "Center of cylindrical boundaries", &cylindricalCenter);
   opts.require ("cylindricalBC", "cylindricalBCl1", "Length of first cylinder",
                 &cylindricalBCl1);
   opts.range("cylindricalBCl1", POSITIVE);
   opts.units("cylindricalBCl1", N_ANGSTROM);
   opts.optional ("cylindricalBCl1", "cylindricalBCl2", "Length of second cylinder",
                  &cylindricalBCl2);
   opts.range ("cylindricalBCl2", POSITIVE);
   opts.units ("cylindricalBCl2", N_ANGSTROM);
// end  additions

   opts.optional("cylindricalBCr1", "cylindricalBCr2", "Radius for second cylinder "
                 "potential", &cylindricalBCr2);
   opts.range("cylindricalBCr2", POSITIVE);
   opts.units("cylindricalBCr2", N_ANGSTROM);
   opts.require("cylindricalBCr2", "cylindricalBCk2", "Force constant for second "
                "cylinder potential (+ is an inward force, - outward)",
                &cylindricalBCk2);
   opts.units("cylindricalBCk2", N_KCAL);
   opts.optional("cylindricalBCr2", "cylindricalBCexp2", "Exponent for second "
                "cylinder potential", &cylindricalBCexp2, 2);
   opts.range("cylindricalBCexp2", POSITIVE);

   ///////////////  Electric field options
   opts.optionalB("main", "eFieldOn", "Should an electric field be applied",
                 &eFieldOn, FALSE);
   opts.optionalB("eFieldOn", "eFieldNormalized", "Is eField vector scaled by cell basis vectors?",
                 &eFieldNormalized, FALSE);
   opts.require("eFieldOn", "eField", "Electric field vector", &eField);
   opts.optional("eFieldOn", "eFieldFreq", "Electric field frequency", &eFieldFreq);
   opts.optional("eFieldOn", "eFieldPhase", "Electric field phase", &eFieldPhase);

      ///////////////  Stir options
   opts.optionalB("main", "stirOn", "Should stirring torque be applied",
                 &stirOn, FALSE);
   opts.optional("stirOn", "stirFilename", "PDB file with flags for "
     "stirred atoms (default is the PDB input file)",
		 PARSE_STRING);
   opts.optional("stirOn", "stirredAtomsCol", "Column in the stirredAtomsFile "
		 "containing the flags (nonzero means fixed);\n"
		 "default is 'O'", PARSE_STRING);
   opts.require("stirOn", "stirStartingTheta", "Stir starting theta offset", &stirStartingTheta);
   opts.require("stirOn", "stirK", "Stir force harmonic spring constant", &stirK);
   //should make this optional, compute from firsttimestep * stirVel
   opts.require("stirOn", "stirVel", "Stir angular velocity (deg/timestep)", &stirVel);
   opts.require("stirOn", "stirAxis", "Stir axis (direction vector)", &stirAxis);
   opts.require("stirOn", "stirPivot", "Stir pivot point (coordinate)", &stirPivot);

    //////////  Extra bonds options
   opts.optionalB("main", "extraBonds",
		"Should extra bonded forces be applied",
                 &extraBondsOn, FALSE);
   opts.optional("extraBonds", "extraBondsFile",
		"file with list of extra bonds",
		 PARSE_MULTIPLES);
   opts.optionalB("extraBonds", "extraBondsCosAngles",
		"Should extra angles be cosine-based to match ancient bug",
		&extraBondsCosAngles, TRUE);

}

void SimParameters::config_parser_misc(ParseOptions &opts) {

   ///////////////  Load balance options
   opts.optional("main", "ldBalancer", "Load balancer",
     loadBalancer);
   opts.optional("main", "ldbStrategy", "Load balancing strategy",
     loadStrategy);
   opts.optional("main", "ldbPeriod", "steps between load balancing",
     &ldbPeriod);
   opts.range("ldbPeriod", POSITIVE);
   opts.optional("main", "firstLdbStep", "when to start load balancing",
     &firstLdbStep);
   opts.range("firstLdbStep", POSITIVE);
   opts.optional("main", "lastLdbStep", "when to stop load balancing",
     &lastLdbStep);
   opts.range("lastLdbStep", POSITIVE);
   opts.optional("main", "hybridGroupSize", "Hybrid load balancing group size",
     &hybridGroupSize);
   opts.optional("main", "ldbBackgroundScaling",
     "background load scaling", &ldbBackgroundScaling);
   opts.range("ldbBackgroundScaling", NOT_NEGATIVE);
   opts.optional("main", "ldbPMEBackgroundScaling",
     "PME node background load scaling", &ldbPMEBackgroundScaling);
   opts.range("ldbPMEBackgroundScaling", NOT_NEGATIVE);
   opts.optional("main", "ldbHomeBackgroundScaling",
     "home node background load scaling", &ldbHomeBackgroundScaling);
   opts.range("ldbHomeBackgroundScaling", NOT_NEGATIVE);
   opts.optional("main", "ldbRelativeGrainsize",
     "fraction of average load per compute", &ldbRelativeGrainsize, 0.);
   opts.range("ldbRelativeGrainsize", NOT_NEGATIVE);

   opts.optional("main", "traceStartStep", "when to start tracing", &traceStartStep);
   opts.range("traceStartStep", POSITIVE);
   opts.optional("main", "numTraceSteps", "the number of timesteps to be traced", &numTraceSteps);
   opts.range("numTraceSteps", POSITIVE);

#ifdef MEASURE_NAMD_WITH_PAPI
   opts.optionalB("main", "papiMeasure", "whether use PAPI to measure performacne", &papiMeasure, FALSE);
   opts.optional("main", "papiMeasureStartStep", "when to measure performacne using PAPI", &papiMeasureStartStep);
   opts.range("papiMeasureStartStep", POSITIVE);
   opts.optional("main", "numPapiMeasureSteps", "the number of timesteps to be measured using PAPI", &numPapiMeasureSteps);
   opts.range("numPapiMeasureSteps", POSITIVE);
#endif

   opts.optionalB("main", "outputMaps", "whether to dump compute map and patch map for analysis just before load balancing", &outputMaps, FALSE);
   opts.optionalB("main", "benchTimestep", "whether to do benchmarking timestep in which case final file output is disabled", &benchTimestep, FALSE);
   opts.optional("main", "useCkLoop", "whether to use CkLoop library to parallelize a loop in a function like OpenMP", &useCkLoop,
    #if CMK_SMP && USE_CKLOOP
     ( CkNumPes() < 2 * CkNumNodes() ? 0 : CKLOOP_CTRL_PME_FORWARDFFT ) );
    #else
     0);
    #endif
   opts.range("useCkLoop", NOT_NEGATIVE);

   opts.optionalB("main", "simulateInitialMapping", "whether to study the initial mapping scheme", &simulateInitialMapping, FALSE);
   opts.optional("main", "simulatedPEs", "the number of PEs to be used for studying initial mapping", &simulatedPEs);
   opts.range("simulatedPEs", POSITIVE);
   opts.optional("main", "simulatedNodeSize", "the node size to be used for studying initial mapping", &simulatedNodeSize);
   opts.range("simulatedNodeSize", POSITIVE);
   opts.optionalB("main", "disableTopology", "ignore torus information during patch placement", &disableTopology, FALSE);
   opts.optionalB("main", "verboseTopology", "print torus information during patch placement", &verboseTopology, FALSE);

   opts.optionalB("main", "ldbUnloadPME", "no load on PME nodes",
     &ldbUnloadPME, FALSE);
   opts.optionalB("main", "ldbUnloadZero", "no load on pe zero",
     &ldbUnloadZero, FALSE);
   opts.optionalB("main", "ldbUnloadOne", "no load on pe one",
     &ldbUnloadOne, FALSE);
   opts.optionalB("main", "ldbUnloadOutputPEs", "no load on output PEs",
     &ldbUnloadOutputPEs, FALSE);
   opts.optionalB("main", "noPatchesOnZero", "no patches on pe zero",
     &noPatchesOnZero, FALSE);
   opts.optionalB("main", "noPatchesOnOutputPEs", "no patches on Output PEs",
     &noPatchesOnOutputPEs, FALSE);
   opts.optionalB("main", "noPatchesOnOne", "no patches on pe one",
     &noPatchesOnOne, FALSE);
   opts.optionalB("main", "useCompressedPsf", "The structure file psf is in the compressed format",
                  &useCompressedPsf, FALSE);
   opts.optionalB("main", "genCompressedPsf", "Generate the compressed version of the psf file",
                  &genCompressedPsf, FALSE);
   opts.optionalB("main", "usePluginIO", "Use the plugin I/O to load the molecule system",
                  &usePluginIO, FALSE);
   opts.optionalB("main", "mallocTest", "test how much memory all PEs can allocate",
                  &mallocTest, FALSE);
   opts.optionalB("main", "printExclusions", "print exclusion lists to stdout",
                  &printExclusions, FALSE);
   opts.optional("main", "proxySendSpanningTree", "using spanning tree to send proxies",
                  &proxySendSpanningTree, -1);
   opts.optional("main", "proxyRecvSpanningTree", "using spanning tree to receive proxies",
                  &proxyRecvSpanningTree, 0);  // default off due to memory leak -1);
   opts.optional("main", "proxyTreeBranchFactor", "the branch factor when building a spanning tree",
                  &proxyTreeBranchFactor, 0);  // actual default in ProxyMgr.C
   opts.optionalB("main", "twoAwayX", "half-size patches in 1st dimension",
     &twoAwayX, -1);
   opts.optionalB("main", "twoAwayY", "half-size patches in 2nd dimension",
     &twoAwayY, -1);
   opts.optionalB("main", "twoAwayZ", "half-size patches in 3rd dimension",
     &twoAwayZ, -1);
   opts.optional("main", "maxPatches", "maximum patch count", &maxPatches, -1);

   /////  Restart timestep option
   opts.optional("main", "firsttimestep", "Timestep to start simulation at",
     &firstTimestep, 0);
   opts.range("firsttimestep", NOT_NEGATIVE);

   /////  Test mode options
   opts.optionalB("main", "test", "Perform self-tests rather than simulation",
		&testOn, FALSE);
   opts.optionalB("main", "commOnly", "Do not evaluate forces or integrate",
		&commOnly, FALSE);

   opts.optionalB("main", "statsOn", "counters in machine layer",
		&statsOn, FALSE);
   ///////////////  hydrogen bond computation options
   opts.optionalB("main", "hbonds", "Use explicit hydrogen bond term",
                 &HydrogenBonds, FALSE);
   opts.optionalB("hbonds","hbAntecedents","Include Antecedent in hbond term",
                 &useAntecedent, TRUE);
   opts.optional("hbonds","hbAAexp","Hbond AA-A-H angle cos exponential",
                 &aaAngleExp, 2);
   opts.optional("hbonds","hbHAexp","Hbond D-H-A angle cos exponential",
                 &haAngleExp, 4);
   opts.optional("hbonds","hbDistAexp","Hbond A-D dist attractive exponential",
                 &distAttExp, 4);
   opts.optional("hbonds","hbDistRexp","Hbond A-D dist repulstive exponential",
                 &distRepExp, 6);
   opts.optional("hbonds","hbCutoffAngle","Hbond D-H-A cutoff angle",
                 &dhaCutoffAngle, 100.0);
   opts.range("hbCutoffAngle", NOT_NEGATIVE);
   opts.optional("hbonds","hbOnAngle","Hbond D-H-A switch function on angle",
                 &dhaOnAngle, 60.0);
   opts.range("hbOnAngle", NOT_NEGATIVE);
   opts.optional("hbonds","hbOffAngle","Hbond D-H-A switch function off angle",
                 &dhaOffAngle, 80.0);
   opts.range("hbOffAngle", NOT_NEGATIVE);
   opts.optional("hbonds","hbCutoffDist","Hbond A-D cutoff distance",
                 &daCutoffDist, 7.5);
   opts.range("hbCutoffDist", POSITIVE);
   opts.units("hbCutoffDist", N_ANGSTROM);
   opts.optional("hbonds","hbOnDist","Hbond A-D switch function on distance",
                 &daOnDist, 5.5);
   opts.range("hbOnDist", POSITIVE);
   opts.units("hbOnDist", N_ANGSTROM);
   opts.optional("hbonds","hbOffDist","Hbond A-D switch function off distance",
                 &daOffDist, 6.5);
   opts.range("hbOffDist", POSITIVE);
   opts.units("hbOffDist", N_ANGSTROM);

   // IMD options
   opts.optionalB("main","IMDon","Connect using IMD?",&IMDon, FALSE);
   opts.require("IMDon","IMDport", "Port to which to bind", &IMDport);
   opts.range("IMDport",POSITIVE);
   opts.require("IMDon","IMDfreq", "Frequency at which to report", &IMDfreq);
   opts.range("IMDfreq",POSITIVE);
   opts.optionalB("IMDon","IMDwait","Pause until IMD connection?",&IMDwait,
     FALSE);
   opts.optionalB("IMDon","IMDignore","Ignore any user input?",&IMDignore,
     FALSE);
   opts.optionalB("IMDon","IMDignoreForces","Ignore forces ONLY?",&IMDignoreForces,
     FALSE);
   // Maximum Partition options
   opts.optional("ldBalancer", "maxSelfPart",
     "maximum number of self partitions in one patch", &maxSelfPart, 20);
   opts.range("maxSelfPart",POSITIVE);
   opts.optional("ldBalancer", "maxPairPart",
     "maximum number of pair partitions in one patch", &maxPairPart, 8);
   opts.range("maxPairPart",POSITIVE);
   opts.optional("ldBalancer", "numAtomsSelf",
		 "maximum number of atoms in one self compute distribution",
		 &numAtomsSelf, 154);
   opts.range("numAtomsSelf",NOT_NEGATIVE);

   opts.optional("ldBalancer", "numAtomsSelf2",
		 "maximum number of atoms in one self compute distribution",
		 &numAtomsSelf2, 154);
   opts.range("numAtomsSelf2",NOT_NEGATIVE);

   opts.optional("ldBalancer", "numAtomsPair",
		 "maximum number of atoms in one pair compute distribution",
		 &numAtomsPair, 318);
   opts.range("numAtomsPair",NOT_NEGATIVE);
   opts.optional("ldBalancer", "numAtomsPair2",
               "maximum number of atoms in one pair compute distribution",
               &numAtomsPair2, 637);
   opts.range("numAtomsPair2",NOT_NEGATIVE);
   opts.optional("main", "minAtomsPerPatch",
               "minimum average atoms per patch",
               &minAtomsPerPatch, 40);
   opts.range("minAtomsPerPatch",NOT_NEGATIVE);

   // Maximum exclusion flags per atom
   opts.optional("main", "maxExclusionFlags",
     "maximum number of exclusion flags per atom", &maxExclusionFlags, 256);
   opts.range("maxExclusionFlags",POSITIVE);

   // Bonded interactions on GPU
   opts.optional("main", "bondedCUDA", "Bitmask for calculating bonded interactions on GPU", &bondedCUDA, 255);

   // Automatically disable individual CUDA kernels that are
   // incompatible with simulation options.
   // Set FALSE to manually control kernel use for development.
   opts.optionalB("main", "useCUDAdisable", "Disable kernels to maintain feature compatibility with CUDA", &useCUDAdisable, TRUE);

   // MIC specific parameters
   opts.optional("main", "mic_unloadMICPEs", "Indicates whether or not the load balancer should unload PEs driving Xeon Phi cards", &mic_unloadMICPEs, 1);
   opts.optional("main", "mic_singleKernel", "Set to non-zero to have all MIC work to be placed in a single kernel", &mic_singleKernel, 1);
   opts.optional("main", "mic_deviceThreshold", "Threshold to use for directing computes to Xeon Phi devices", &mic_deviceThreshold, -1);
   opts.optional("main", "mic_hostSplit", "DMK - reserved", &mic_hostSplit, -1);
   opts.optional("main", "mic_numParts_self_p1", "MIC-Specific NumParts SELF Parameter 1", &mic_numParts_self_p1, -1);
   opts.optional("main", "mic_numParts_pair_p1", "MIC-Specific NumParts PAIR Parameter 1", &mic_numParts_pair_p1, -1);
   opts.optional("main", "mic_numParts_pair_p2", "MIC-Specific NumParts PAIR Parameter 2", &mic_numParts_pair_p2, -1);
   opts.range("mic_unloadMICPEs", NOT_NEGATIVE);
   opts.range("mic_singleKernel", NOT_NEGATIVE);
}

void SimParameters::readExtendedSystem(const char *filename, Lattice *latptr) {

     if ( ! latptr ) {
       iout << iINFO << "EXTENDED SYSTEM FILE   " << filename << "\n" << endi;
     }

     ifstream xscFile(filename);
     if ( ! xscFile ) NAMD_die("Unable to open extended system file.\n");

     char labels[1024];
     do {
       if ( ! xscFile ) NAMD_die("Error reading extended system file.\n");
       xscFile.getline(labels,1023);
     } while ( strncmp(labels,"#$LABELS ",9) );

     int a_x, a_y, a_z, b_x, b_y, b_z, c_x, c_y, c_z;
     a_x = a_y = a_z = b_x = b_y = b_z = c_x = c_y = c_z = -1;
     int o_x, o_y, o_z, s_u, s_v, s_w, s_x, s_y, s_z;
     o_x = o_y = o_z = s_u = s_v = s_w = s_x = s_y = s_z = -1;

     int pos = 0;
     char *l_i = labels + 8;
     while ( *l_i ) {
       if ( *l_i == ' ' ) { ++l_i; continue; }
       char *l_i2;
       for ( l_i2 = l_i; *l_i2 && *l_i2 != ' '; ++l_i2 );
       if ( (l_i2 - l_i) == 3 && (l_i[1] == '_') ) {
	 if (l_i[0] == 'a' && l_i[2] == 'x') a_x = pos;
	 if (l_i[0] == 'a' && l_i[2] == 'y') a_y = pos;
	 if (l_i[0] == 'a' && l_i[2] == 'z') a_z = pos;
	 if (l_i[0] == 'b' && l_i[2] == 'x') b_x = pos;
	 if (l_i[0] == 'b' && l_i[2] == 'y') b_y = pos;
	 if (l_i[0] == 'b' && l_i[2] == 'z') b_z = pos;
	 if (l_i[0] == 'c' && l_i[2] == 'x') c_x = pos;
	 if (l_i[0] == 'c' && l_i[2] == 'y') c_y = pos;
	 if (l_i[0] == 'c' && l_i[2] == 'z') c_z = pos;
	 if (l_i[0] == 'o' && l_i[2] == 'x') o_x = pos;
	 if (l_i[0] == 'o' && l_i[2] == 'y') o_y = pos;
	 if (l_i[0] == 'o' && l_i[2] == 'z') o_z = pos;
	 if (l_i[0] == 's' && l_i[2] == 'u') s_u = pos;
	 if (l_i[0] == 's' && l_i[2] == 'v') s_v = pos;
	 if (l_i[0] == 's' && l_i[2] == 'w') s_w = pos;
	 if (l_i[0] == 's' && l_i[2] == 'x') s_x = pos;
	 if (l_i[0] == 's' && l_i[2] == 'y') s_y = pos;
	 if (l_i[0] == 's' && l_i[2] == 'z') s_z = pos;
       }
       ++pos;
       l_i = l_i2;
     }
     int numpos = pos;

     for ( pos = 0; pos < numpos; ++pos ) {
       double tmp;
       xscFile >> tmp;
       if ( ! xscFile ) NAMD_die("Error reading extended system file.\n");
       if ( pos == a_x ) cellBasisVector1.x = tmp;
       if ( pos == a_y ) cellBasisVector1.y = tmp;
       if ( pos == a_z ) cellBasisVector1.z = tmp;
       if ( pos == b_x ) cellBasisVector2.x = tmp;
       if ( pos == b_y ) cellBasisVector2.y = tmp;
       if ( pos == b_z ) cellBasisVector2.z = tmp;
       if ( pos == c_x ) cellBasisVector3.x = tmp;
       if ( pos == c_y ) cellBasisVector3.y = tmp;
       if ( pos == c_z ) cellBasisVector3.z = tmp;
       if ( pos == o_x ) cellOrigin.x = tmp;
       if ( pos == o_y ) cellOrigin.y = tmp;
       if ( pos == o_z ) cellOrigin.z = tmp;
       if ( pos == s_u ) strainRate2.x = tmp;
       if ( pos == s_v ) strainRate2.y = tmp;
       if ( pos == s_w ) strainRate2.z = tmp;
       if ( pos == s_x ) strainRate.x = tmp;
       if ( pos == s_y ) strainRate.y = tmp;
       if ( pos == s_z ) strainRate.z = tmp;
     }

   if ( latptr ) {
     Lattice test;
     test.set(cellBasisVector1,cellBasisVector2,cellBasisVector3,cellOrigin);

     if ( test.a_p() && ! lattice.a_p() ) {
       NAMD_die("cellBasisVector1 added during atom reinitialization");
     }
     if ( lattice.a_p() && ! test.a_p() ) {
       NAMD_die("cellBasisVector1 dropped during atom reinitialization");
     }
     if ( test.b_p() && ! lattice.b_p() ) {
       NAMD_die("cellBasisVector2 added during atom reinitialization");
     }
     if ( lattice.b_p() && ! test.b_p() ) {
       NAMD_die("cellBasisVector2 dropped during atom reinitialization");
     }
     if ( test.c_p() && ! lattice.c_p() ) {
       NAMD_die("cellBasisVector3 added during atom reinitialization");
     }
     if ( lattice.c_p() && ! test.c_p() ) {
       NAMD_die("cellBasisVector3 dropped during atom reinitialization");
     }

     latptr->set(cellBasisVector1,cellBasisVector2,cellBasisVector3,cellOrigin);
   }

}

#ifdef MEM_OPT_VERSION
//This global var is defined in mainfunc.C
extern char *gWorkDir;
#endif

void SimParameters::check_config(ParseOptions &opts, ConfigList *config, char *&cwd) {

   int len;    //  String length
   StringList *current; //  Pointer to config option list

#ifdef MEM_OPT_VERSION
   char *namdWorkDir = NULL;
#endif

  if ( opts.defined("obsolete") ) {
    iout << iWARN <<
      "\"obsolete\" defined, silently ignoring obsolete options\n" << endi;
  }

   //  Take care of cwd processing
   if (opts.defined("cwd"))
   {
    //  First allocate and get the cwd value
    current = config->find("cwd");

    len = strlen(current->data);

    if ( CHDIR(current->data) )
    {
      NAMD_die("chdir() to given cwd failed!");
    } else {
      iout << iINFO << "Changed directory to " << current->data << "\n" << endi;
    }

    if (current->data[len-1] != PATHSEP)
      len++;

    cwd = new char[len+1];

    strcpy(cwd, current->data);

    if (current->data[strlen(current->data)-1] != PATHSEP)
      strcat(cwd, PATHSEPSTR);
   }

#ifdef MEM_OPT_VERSION
   if(cwd!=NULL)namdWorkDir = cwd;
   else namdWorkDir = gWorkDir;
   int dirlen = strlen(namdWorkDir);
   //only support the path representation on UNIX-like platforms
   char *tmpDir;
   if(namdWorkDir[dirlen-1]=='/'){
     tmpDir = new char[dirlen+1];
     tmpDir[dirlen] = 0;
   }else{
     tmpDir = new char[dirlen+2];
     tmpDir[dirlen]='/';
     tmpDir[dirlen+1]=0;
   }
   memcpy(tmpDir, namdWorkDir, dirlen);
   namdWorkDir = tmpDir;
 //finished recording the per atom files, free the space for gWorkDir
   delete [] gWorkDir;
#endif


   // Don't try to specify coordinates with pluginIO
   if ( usePluginIO && opts.defined("coordinates") ) {
     NAMD_die("Separate coordinates file not allowed with plugin IO, coordinates will be taken from structure file.");
   }

   // If it's not AMBER||GROMACS, then "coordinates", "structure"
   // and "parameters" must be specified.
   if (!amberOn && !gromacsOn) {
#ifndef MEM_OPT_VERSION
     if (useCompressedPsf)
       NAMD_die("useCompressedPsf requires memory-optimized build!");
     if (!usePluginIO && !genCompressedPsf && !opts.defined("coordinates"))
       NAMD_die("coordinates not found in the configuration file!");
#else
     if(!usePluginIO && !opts.defined("bincoordinates")) {
       NAMD_die("bincoordinates not found in the configuration file for the memory optimized version!");
     }
     if(!usePluginIO && opts.defined("coordinates")) {
       NAMD_die("coordinates not allowed in the configuration file for the memory optimized version!");
     }
#endif
     if (!opts.defined("structure"))
       NAMD_die("structure not found in the configuration file!");
     if (!opts.defined("parameters"))
       NAMD_die("parameters not found in the configuration file!");
   }

   // In any case, there should be either "coordinates" or
   // "ambercoor", but not both
   if (opts.defined("coordinates") && opts.defined("ambercoor"))
     NAMD_die("Cannot specify both coordinates and ambercoor!");
#ifndef MEM_OPT_VERSION
   if (!genCompressedPsf && !opts.defined("coordinates") && !opts.defined("ambercoor")
       && !opts.defined("grocoorfile") && !usePluginIO)
     NAMD_die("Coordinate file not found!");
#endif

   //  Make sure that both a temperature and a velocity PDB were
   //  specified
   if (opts.defined("temperature") &&
       (opts.defined("velocities") || opts.defined("binvelocities")) )
   {
      NAMD_die("Cannot specify both an initial temperature and a velocity file");
   }

#ifdef MEM_OPT_VERSION
//record the absolute file name for binAtomFile, binCoorFile and binVelFile etc.
   binAtomFile = NULL;
   binCoorFile = NULL;
   binVelFile = NULL;
   binRefFile = NULL;

   char *curfile = NULL;
   dirlen = strlen(namdWorkDir);
   current = config->find("structure");;
   curfile = current->data;
   int filelen = strlen(curfile);
   if(*curfile == '/' || *curfile=='~') {
     //check whether it is an absolute path
     //WARNING: Only works on Unix-like platforms!
     //Needs to fix on Windows platform.
     //-Chao Mei
     //adding 5 because of ".bin"+"\0"
     binAtomFile = new char[filelen+5];
     memcpy(binAtomFile, curfile, filelen);
     memcpy(binAtomFile+filelen, ".bin", 4);
     binAtomFile[filelen+4] = 0;
   }else{
     binAtomFile = new char[dirlen+filelen+5];
     memcpy(binAtomFile, namdWorkDir, dirlen);
     memcpy(binAtomFile+dirlen, curfile, filelen);
     memcpy(binAtomFile+dirlen+filelen, ".bin", 4);
     binAtomFile[dirlen+filelen+4] = 0;
   }

   current = config->find("bincoordinates");
   curfile = current->data;
   filelen = strlen(curfile);
   if(*curfile == '/' || *curfile=='~') {
     binCoorFile = new char[filelen+1];
     memcpy(binCoorFile, curfile, filelen);
     binCoorFile[filelen] = 0;
   }else{
     binCoorFile = new char[dirlen+filelen+1];
     memcpy(binCoorFile, namdWorkDir, dirlen);
     memcpy(binCoorFile+dirlen, curfile, filelen);
     binCoorFile[dirlen+filelen] = 0;
   }

   if(opts.defined("binvelocities")){
     current = config->find("binvelocities");
     curfile = current->data;
     filelen = strlen(curfile);
     if(*curfile == '/' || *curfile=='~') {
       binVelFile = new char[filelen+1];
       memcpy(binVelFile, curfile, filelen);
       binVelFile[filelen] = 0;
     }else{
       binVelFile = new char[dirlen+filelen+1];
       memcpy(binVelFile, namdWorkDir, dirlen);
       memcpy(binVelFile+dirlen, curfile, filelen);
       binVelFile[dirlen+filelen] = 0;
     }
   }

   if(opts.defined("binrefcoords")){
     current = config->find("binrefcoords");
     curfile = current->data;
     filelen = strlen(curfile);
     if(*curfile == '/' || *curfile=='~') {
       binRefFile = new char[filelen+1];
       memcpy(binRefFile, curfile, filelen);
       binRefFile[filelen] = 0;
     }else{
       binRefFile = new char[dirlen+filelen+1];
       memcpy(binRefFile, namdWorkDir, dirlen);
       memcpy(binRefFile+dirlen, curfile, filelen);
       binRefFile[dirlen+filelen] = 0;
     }
   }

   //deal with output file name to make it absolute path for parallel output
   if(outputFilename[0] != '/' && outputFilename[0]!='~') {
     filelen = strlen(outputFilename);
     char *tmpout = new char[filelen];
     memcpy(tmpout, outputFilename, filelen);
     CmiAssert(filelen+dirlen <= 120); //leave 8 chars for file suffix
     memcpy(outputFilename, namdWorkDir, dirlen);
     memcpy(outputFilename+dirlen, tmpout, filelen);
     outputFilename[filelen+dirlen] = 0;
     delete [] tmpout;
   }

   if ( dcdFrequency && opts.defined("dcdfile") &&
        dcdFilename[0] != '/' && dcdFilename[0]!='~' ) {
     filelen = strlen(dcdFilename);
     char *tmpout = new char[filelen];
     memcpy(tmpout, dcdFilename, filelen);
     CmiAssert(filelen+dirlen <= 120); //leave 8 chars for file suffix
     memcpy(dcdFilename, namdWorkDir, dirlen);
     memcpy(dcdFilename+dirlen, tmpout, filelen);
     dcdFilename[filelen+dirlen] = 0;
     delete [] tmpout;
   }

   if ( velDcdFrequency && opts.defined("veldcdfile") &&
        velDcdFilename[0] != '/' && velDcdFilename[0]!='~' ) {
     filelen = strlen(velDcdFilename);
     char *tmpout = new char[filelen];
     memcpy(tmpout, velDcdFilename, filelen);
     CmiAssert(filelen+dirlen <= 120); //leave 8 chars for file suffix
     memcpy(velDcdFilename, namdWorkDir, dirlen);
     memcpy(velDcdFilename+dirlen, tmpout, filelen);
     velDcdFilename[filelen+dirlen] = 0;
     delete [] tmpout;
   }

   if ( forceDcdFrequency && opts.defined("forcedcdfile") &&
        forceDcdFilename[0] != '/' && forceDcdFilename[0]!='~' ) {
     filelen = strlen(forceDcdFilename);
     char *tmpout = new char[filelen];
     memcpy(tmpout, forceDcdFilename, filelen);
     CmiAssert(filelen+dirlen <= 120); //leave 8 chars for file suffix
     memcpy(forceDcdFilename, namdWorkDir, dirlen);
     memcpy(forceDcdFilename+dirlen, tmpout, filelen);
     forceDcdFilename[filelen+dirlen] = 0;
     delete [] tmpout;
   }

   if ( restartFrequency && opts.defined("restartname") &&
        restartFilename[0] != '/' && restartFilename[0]!='~' ) {
     filelen = strlen(restartFilename);
     char *tmpout = new char[filelen];
     memcpy(tmpout, restartFilename, filelen);
     CmiAssert(filelen+dirlen <= 120); //leave 8 chars for file suffix
     memcpy(restartFilename, namdWorkDir, dirlen);
     memcpy(restartFilename+dirlen, tmpout, filelen);
     restartFilename[filelen+dirlen] = 0;
     delete [] tmpout;
   }

   delete [] namdWorkDir;

   if (opts.defined("numinputprocs")) {
     if(numinputprocs > CkNumPes()) {
       iout << iWARN << "The number of input processors exceeds the total number of processors. Resetting to half of the number of total processors.\n" << endi;
       numinputprocs = (CkNumPes()>>1)+(CkNumPes()&1);
     }
   }

   if (opts.defined("numoutputprocs")) {
     if(numoutputprocs > CkNumPes()) {
       iout << iWARN << "The number of output processors exceeds the total number of processors. Resetting to half of the number of total processors.\n" << endi;
       numoutputprocs = (CkNumPes()>>1)+(CkNumPes()&1);
     }
   }

#ifndef OUTPUT_SINGLE_FILE
#error OUTPUT_SINGLE_FILE not defined!
#endif

   #if !OUTPUT_SINGLE_FILE
   //create directories for multi-file output scheme
   create_output_directories("coor");
   create_output_directories("vel");
   if(dcdFrequency) {
	   create_output_directories("dcd");
	   if(opts.defined("dcdfile")){
		   iout << iWARN << "The dcd file output has been changed to directory: " << outputFilename << ".\n" << endi;
	   }
   }
   if (velDcdFrequency) {
	   create_output_directories("veldcd");
	   if(opts.defined("veldcdfile")){
		   iout << iWARN << "The veldcd file output has been changed to directory: " << outputFilename << ".\n" << endi;
	   }
   }
   if (forceDcdFrequency) {
	   create_output_directories("forcedcd");
	   if(opts.defined("forcedcdfile")){
		   iout << iWARN << "The forcedcd file output has been changed to directory: " << outputFilename << ".\n" << endi;
	   }
   }
   #endif
#endif

   if (! opts.defined("auxFile")) {
     strcpy(auxFilename,outputFilename);
     strcat(auxFilename,".aux");
   }

   //  Check for frequencies
   if (dcdFrequency) {
     if (! opts.defined("dcdfile")) {
       strcpy(dcdFilename,outputFilename);
       strcat(dcdFilename,".dcd");
     }
   } else {
     dcdFilename[0] = STRINGNULL;
   }

   if (velDcdFrequency) {
     if (! opts.defined("veldcdfile")) {
       strcpy(velDcdFilename,outputFilename);
       strcat(velDcdFilename,".veldcd");
     }
   } else {
     velDcdFilename[0] = STRINGNULL;
   }

   if (forceDcdFrequency) {
     if (! opts.defined("forcedcdfile")) {
       strcpy(forceDcdFilename,outputFilename);
       strcat(forceDcdFilename,".forcedcd");
     }
   } else {
     forceDcdFilename[0] = STRINGNULL;
   }

   if (xstFrequency) {
     if (! opts.defined("xstfile")) {
       strcpy(xstFilename,outputFilename);
       strcat(xstFilename,".xst");
     }
   } else {
     xstFilename[0] = STRINGNULL;
   }

   if (restartFrequency) {
     if (! opts.defined("restartname")) {
       strcpy(restartFilename,outputFilename);
       if ( ! restartSave ) strcat(restartFilename,".restart");
     }
   } else {
     restartFilename[0] = STRINGNULL;
     restartSave = FALSE;
     binaryRestart = FALSE;
   }

   if (storeComputeMap || loadComputeMap) {
     if (! opts.defined("computeMapFile")) {
       strcpy(computeMapFilename,"computeMapFile");
       strcat(computeMapFilename,".txt");
     }
   }


   if (!amberOn)
   { //****** BEGIN CHARMM/XPLOR type changes
     //// set default
     if (!paraTypeXplorOn && !paraTypeCharmmOn)
     {
       paraTypeXplorOn = TRUE;
     }
     //// make sure that there is just one type of input parameters specified
     if (paraTypeXplorOn && paraTypeCharmmOn)
     {
       NAMD_die("Please specify either XPLOR or CHARMM format for parameters!");
     }
     //****** END CHARMM/XPLOR type changes
   }


   //  If minimization isn't on, must have a temp or velocity
   if (!(minimizeOn||minimizeCGOn) && !opts.defined("temperature") &&
       !opts.defined("velocities") && !opts.defined("binvelocities") )
   {
      NAMD_die("Must have either an initial temperature or a velocity file");
   }

   if (minimizeOn||minimizeCGOn) { initialTemp = 0.0; }
   if (opts.defined("velocities") || opts.defined("binvelocities") )
   {
  initialTemp = -1.0;
   }

   ///// periodic cell parameters

   if ( opts.defined("extendedSystem") ) readExtendedSystem(config->find("extendedSystem")->data);

#ifdef MEM_OPT_VERSION
   if ( LJcorrection ) {
      NAMD_die("LJ tail corrections not yet available for memory optimized builds");
   }
#endif

   if ( LJcorrection && ! cellBasisVector3.length2() ) {
     NAMD_die("Can't use LJ tail corrections without periodic boundary conditions!");
   }

   if ( cellBasisVector3.length2() && ! cellBasisVector2.length2() ) {
     NAMD_die("Used cellBasisVector3 without cellBasisVector2!");
   }

   if ( cellBasisVector2.length2() && ! cellBasisVector1.length2() ) {
     NAMD_die("Used cellBasisVector2 without cellBasisVector1!");
   }

   if ( cellOrigin.length2() && ! cellBasisVector1.length2() ) {
     NAMD_die("Used cellOrigin without cellBasisVector1!");
   }

   lattice.set(cellBasisVector1,cellBasisVector2,cellBasisVector3,cellOrigin);

   if (! opts.defined("DCDunitcell")) {
      dcdUnitCell = lattice.a_p() && lattice.b_p() && lattice.c_p();
   }

   char s[129];

   ///// cylindricalBC stuff
   if ( ! opts.defined("cylindricalBCAxis") )
   {
      cylindricalBCAxis = 'x';
   }
   else
   {
     opts.get("cylindricalBCAxis", s);

     if (!strcasecmp(s, "x"))
     {
      cylindricalBCAxis = 'x';
     }
     else if (!strcasecmp(s, "y"))
     {
      cylindricalBCAxis = 'y';
     }
     else if (!strcasecmp(s, "z"))
     {
      cylindricalBCAxis = 'z';
     }
   else
     {
      char err_msg[128];

      sprintf(err_msg, "Illegal value '%s' for 'cylindricalBCAxis' in configuration file", s);
      NAMD_die(err_msg);
     }
   }

   if (!opts.defined("splitPatch"))
   {
     splitPatch = SPLIT_PATCH_HYDROGEN;
   }
   else
   {
     opts.get("splitPatch", s);
     if (!strcasecmp(s, "position"))
       splitPatch = SPLIT_PATCH_POSITION;
     else if (!strcasecmp(s,"hydrogen"))
       splitPatch = SPLIT_PATCH_HYDROGEN;
     else
     {
       char err_msg[129];
       sprintf(err_msg,
          "Illegal value '%s' for 'splitPatch' in configuration file",
       s);
       NAMD_die(err_msg);
     }
   }

   ///// exclude stuff
   opts.get("exclude", s);

   if (!strcasecmp(s, "none"))
   {
      exclude = NONE;
      splitPatch = SPLIT_PATCH_POSITION;
   }
   else if (!strcasecmp(s, "1-2"))
   {
      exclude = ONETWO;
      splitPatch = SPLIT_PATCH_POSITION;
   }
   else if (!strcasecmp(s, "1-3"))
   {
      exclude = ONETHREE;
   }
   else if (!strcasecmp(s, "1-4"))
   {
      exclude = ONEFOUR;
   }
   else if (!strcasecmp(s, "scaled1-4"))
   {
      exclude = SCALED14;
   }
   else
   {
      char err_msg[128];

      sprintf(err_msg, "Illegal value '%s' for 'exclude' in configuration file",
   s);
      NAMD_die(err_msg);
   }

   if (scale14 != 1.0 && exclude != SCALED14)
   {
      iout << iWARN << "Exclude is not scaled1-4; 1-4scaling ignored.\n" << endi;
   }

   // water model stuff
   if (!opts.defined("waterModel")) {
     watmodel = WAT_TIP3;
   } else {
     opts.get("waterModel", s);
     if (!strncasecmp(s, "tip4", 4)) {
       iout << iINFO << "Using TIP4P water model.\n" << endi;
       watmodel = WAT_TIP4;
     } else if (!strncasecmp(s, "tip3", 4)) {
       iout << iINFO << "Using TIP3P water model.\n" << endi;
       watmodel = WAT_TIP3;
     } else if (!strncasecmp(s, "swm4", 4)) {
       iout << iINFO << "Using SWM4-DP water model.\n" << endi;
       watmodel = WAT_SWM4;
     } else {
       char err_msg[128];
       sprintf(err_msg,
           "Illegal value %s for 'waterModel' in configuration file", s);
       NAMD_die(err_msg);
     }
   }
   if (watmodel == WAT_SWM4 && !drudeOn) {
     NAMD_die("Must have 'drudeOn' enabled to use SWM4-DP water model.");
   }
   if (drudeOn && watmodel != WAT_SWM4) {
     watmodel = WAT_SWM4;
     iout << iWARN
       << "Setting water model to 'swm4' (SWM4-DP) for Drude polarization.\n"
       << endi;
   }

   // Drude water model uses "lonepairs"
   if (watmodel == WAT_SWM4) {
     lonepairs = TRUE;
   }

   // Added by JLai -- 8.2.11 -- Checks if Go method is defined
   if (goForcesOn) {
     iout << iINFO << "Go forces are on\n" << endi;
   // Added by JLai -- 6.3.11 -- Checks if Go method is defined
     int * gomethod = &goMethod;
     if (!opts.defined("GoMethod")) {
       *gomethod = 0;
       //     printf("GO METHOD IS NOT DEFINED SO WE'LL SET IT TO SOME WEIRD VALUE\n");
     } else {
       opts.get("GoMethod",s);
       // printf("GO METHOD IS DEFINED SO WE'LL PRINT IT OUT: %s\n",s);
       *gomethod = atoi(s);
     }
     if (!strcasecmp(s, "matrix")) {
       goMethod = 1;
       //GoMethod = GO_MATRIX;
     } else if (!strcasecmp(s, "faster")) {
       goMethod = 2;
       //GoMethod = GO_FASTER;
     } else if (!strcasecmp(s, "lowmem")) {
       goMethod = 3;
       //GoMethod = GO_LOWMEM;
     }
     else {
       char err_msg[129];
       sprintf(err_msg,
	       "Illegal value '%s' for 'GoMethod' in configuration file",
	       s);
       NAMD_die(err_msg);
     }
   }  // End of NAMD code to check goMethod
   // End of Port -- JL

   //  Get multiple timestep integration scheme
   if (!opts.defined("MTSAlgorithm"))
   {
  MTSAlgorithm = VERLETI;
   }
   else
   {
  opts.get("MTSAlgorithm", s);

  if (!strcasecmp(s, "naive"))
  {
    MTSAlgorithm = NAIVE;
  }
  else if (!strcasecmp(s, "constant"))
  {
    MTSAlgorithm = NAIVE;
  }
  else if (!strcasecmp(s, "impulse"))
  {
    MTSAlgorithm = VERLETI;
  }
  else if (!strcasecmp(s, "verleti"))
  {
    MTSAlgorithm = VERLETI;
  }
  else
  {
    char err_msg[129];

    sprintf(err_msg,
       "Illegal value '%s' for 'MTSAlgorithm' in configuration file",
       s);
    NAMD_die(err_msg);
  }
   }

   //  Get the long range force splitting specification
   if (!opts.defined("longSplitting"))
   {
  longSplitting = C1;
   }
   else
   {
  opts.get("longSplitting", s);
  if (!strcasecmp(s, "sharp"))
    longSplitting = SHARP;
  else if (!strcasecmp(s, "xplor"))
    longSplitting = XPLOR;
  else if (!strcasecmp(s, "c1"))
    longSplitting = C1;
  else if (!strcasecmp(s, "c2"))
    longSplitting = C2;
  else
  {
    char err_msg[129];

    sprintf(err_msg,
       "Illegal value '%s' for 'longSplitting' in configuration file",
       s);
    NAMD_die(err_msg);
  }
   }

   // take care of rigid bond options
   if (!opts.defined("rigidBonds"))
   {
      rigidBonds = RIGID_NONE;
   }
   else
   {
      opts.get("rigidBonds", s);
      if (!strcasecmp(s, "all"))
      {
          rigidBonds = RIGID_ALL;
      }
      else if (!strcasecmp(s, "water"))
      {
           rigidBonds = RIGID_WATER;
      }
      else if (!strcasecmp(s, "none"))
      {
           rigidBonds = RIGID_NONE;
      }
      else
      {
        char err_msg[256];
        sprintf(err_msg,
          "Illegal value '%s' for 'rigidBonds' in configuration file", s);
        NAMD_die(err_msg);
      }
   }

   // TIP4P and SWM4-DP water models require rigid water
   if ((watmodel == WAT_TIP4 || watmodel == WAT_SWM4)
       && rigidBonds == RIGID_NONE) {
     char err_msg[256];
     sprintf(err_msg,
         "Water model %s requires rigidBonds set to \"all\" or \"water\"",
         (watmodel == WAT_TIP4 ? "TIP4P" : "SWM4-DP"));
     NAMD_die(err_msg);
   }

   //  Take care of switching stuff
   if (switchingActive)
   {

     if (!opts.defined("switchDist")) {
       NAMD_die("switchDist must be defined when switching is enabled");
     }

     if ( (switchingDist>cutoff) || (switchingDist<0) )
     {
       char err_msg[129];

       sprintf(err_msg,
         "switchDist muct be between 0 and cutoff, which is %f", cutoff);
       NAMD_die(err_msg);
     }

   }

   if ( martiniSwitching )
   {
     if ( ! switchingActive )
     {
       NAMD_die("martiniSwitching requires switching");
     }
     if ( vdwForceSwitching )
     {
       NAMD_die("martiniSwitching and vdwForceSwitching are exclusive to one another. Select only one.");
     }
     if ( dielectric != 15.0 && ! martiniDielAllow )
     {
       iout << iWARN << "USE DIELECTRIC OF 15.0 WITH MARTINI.\n";
       iout << iWARN << "SETTING dielectric 15.0\n";
       iout << iWARN << "FOR NON-STANDARD DIELECTRIC WITH MARTINI, SET: martiniDielAllow on\n";
       dielectric = 15.0;
     }
     if ( ! cosAngles )
     {
       iout << iWARN << "USE COSINE BASED ANGLES WITH MARTINI.\n";
       iout << iWARN << "SETTING cosAngles on\n";
       cosAngles = TRUE;
     }
     if ( PMEOn )
     {
       NAMD_die("Do not use Particle Mesh Ewald with Martini.  Set: PME off");
     }
     if ( MSMOn )
     {
       NAMD_die("Do not use Multilevel Summation Method with Martini.  Set: MSM off");
     }
     if ( FMMOn )
     {
       NAMD_die("Do not use Fast Multipole Method with Martini.  Set: FMM off");
     }

   }


   if (!opts.defined("pairlistDist"))
   {
  pairlistDist = cutoff;
   }
   else if (pairlistDist < cutoff)
   {
  NAMD_die("pairlistDist must be >= cutoff distance");
   }

   patchDimension = pairlistDist;

   if ( splitPatch == SPLIT_PATCH_HYDROGEN ) {
     patchDimension += hgroupCutoff;
   }

   BigReal defaultMargin = 0.0;
   if (berendsenPressureOn || langevinPistonOn) {
      defaultMargin = ( useFlexibleCell ? 0.06 : 0.03 ) * patchDimension;
   }
   if ( margin == XXXBIGREAL ) {
     margin = defaultMargin;
   }
   if ( defaultMargin != 0.0 && margin == 0.0 ) {
     margin = defaultMargin;
     iout << iWARN << "ALWAYS USE NON-ZERO MARGIN WITH CONSTANT PRESSURE!\n";
     iout << iWARN << "CHANGING MARGIN FROM 0 to " << margin << "\n" << endi;
   }

   patchDimension += margin;

    //ensure patch can handle alpha_cutoff for gbis
    if (GBISOn) {
      //Check compatibility
      if (fullDirectOn) {
        NAMD_die("GBIS not compatible with FullDirect");
      }
      if (PMEOn) {
        NAMD_die("GBIS not compatible with PME");
      }
      if (MSMOn) {
        NAMD_die("GBIS not compatible with MSM");
      }
      if (FMMOn) {
        NAMD_die("GBIS not compatible with FMM");
      }
      if (alchOn) {
        NAMD_die("GBIS not compatible with Alchemical Transformations");
      }
      if (lesOn) {
        NAMD_die("GBIS not compatible with Locally Enhanced Sampling");
      }
      if (FMAOn) {
        NAMD_die("GBIS not compatible with FMA");
      }
      if (drudeOn) {
        NAMD_die("GBIS not compatible with Drude Polarization");
      }

      if (alpha_cutoff > patchDimension) {
        patchDimension = alpha_cutoff;
      }
      //calculate kappa
      BigReal tmp = (initialTemp > 0) ? initialTemp : 300;
      kappa = 50.29216*sqrt(ion_concentration/solvent_dielectric/tmp);
      /*magic number = 1/sqrt(eps0*kB/(2*nA*e^2*1000))*/
    } // GBISOn

    if (LCPOOn) {
#ifdef MEM_OPT_VERSION
      NAMD_die("SASA not yet available for memory optimized builds");
#endif
      if ( lattice.volume() > 0 ) {
        NAMD_die("SASA does not yet support periodic boundary conditions.");
      }
      //LCPO requires patches to be at least 16.2Ang in each dimension
      // twoAway[XYZ} is ignored for now
    }

   //  Turn on global integration if not explicitly specified

   if ( dihedralOn ) globalOn = TRUE;

#ifdef NAMD_CUDA
   if (loweAndersenOn) {
       NAMD_die("Lowe-Andersen dynamics not compatible with CUDA at this time");
   }
#endif
   // END LA

   // BEGIN LA
   if (loweAndersenOn && (langevinOn || tCoupleOn))
   {
      NAMD_die("Lowe-Andersen dynamics, Langevin dynamics and temperature coupling are mutually exclusive dynamics modes");
   }
   // END LA

   if (tCoupleOn && opts.defined("rescaleFreq") )
   {
      NAMD_die("Temperature coupling and temperature rescaling are mutually exclusive");
   }

   if (globalOn && CkNumPes() > 1)
   {
      NAMD_die("Global integration does not run in parallel (yet).");
   }

   if (COLDOn && langevinOn)
   {
      NAMD_die("COLD and Langevin dynamics are mutually exclusive dynamics modes");
   }
   if (COLDOn && minimizeOn)
   {
      NAMD_die("COLD and minimization are mutually exclusive dynamics modes");
   }
   if (COLDOn && tCoupleOn)
   {
      NAMD_die("COLD and temperature coupling are mutually exclusive dynamics modes");
   }
   if (COLDOn && opts.defined("rescaleFreq"))
   {
      NAMD_die("COLD and velocity rescaling are mutually exclusive dynamics modes");
   }

   if (splitPatch == SPLIT_PATCH_POSITION && mollyOn )
   {
      NAMD_die("splitPatch hydrogen is required for MOLLY");
   }

   if (splitPatch == SPLIT_PATCH_POSITION && rigidBonds != RIGID_NONE)
   {
      NAMD_die("splitPatch hydrogen is required for rigidBonds");
   }

   if (accelMDOn) {
       if(accelMDG){
	   char msg[128];
	   if(accelMDGiE < 1 || accelMDGiE > 2){
	       sprintf(msg, "accelMDGiE was set to %d but it should be 1 or 2", accelMDGiE);
	       NAMD_die(msg);
	   }
           if(accelMDGStatWindow > 0){
               if(accelMDGcMDPrepSteps % accelMDGStatWindow != 0)
                   NAMD_die("'accelMDGcMDPrepSteps' has to be a multiple of 'accelMDGStatWindow'");
               if(accelMDGcMDSteps % accelMDGStatWindow != 0)
                   NAMD_die("'accelMDGcMDSteps' has to be a multiple of 'accelMDGStatWindow'");
               if(accelMDGEquiPrepSteps % accelMDGStatWindow != 0)
                   NAMD_die("'accelMDGEquiPrepSteps' has to be a multiple of 'accelMDGStatWindow'");
               if(accelMDGEquiSteps % accelMDGStatWindow != 0)
                   NAMD_die("'accelMDGEquiSteps' has to be a multiple of 'accelMDGStatWindow'");
           }
	   if(accelMDGRestart && accelMDGcMDSteps == 0)
	       accelMDGcMDPrepSteps = 0;
	   else if(accelMDGcMDSteps - accelMDGcMDPrepSteps < 2)
	       NAMD_die("'accelMDGcMDSteps' should be larger than 'accelMDGcMDPrepSteps'");

	   if(accelMDGEquiSteps == 0)
	       accelMDGEquiPrepSteps = 0;
	   else if(accelMDGresetVaftercmd){
	       if(accelMDGEquiPrepSteps <= 0)
		   NAMD_die("'accelMDGEquiPrepSteps' should be non-zero");
	       if(accelMDGEquiSteps - accelMDGEquiPrepSteps < 1)
		   NAMD_die("'accelMDGEquiSteps' should be larger than 'accelMDGEquiPrepSteps'");
	   }

	   //warn user that accelMD params will be ignored
	   if(opts.defined("accelMDE"))
	       iout << iWARN << "accelMDE will be ignored with accelMDG on.\n" << endi;
	   if(opts.defined("accelMDalpha"))
	       iout << iWARN << "accelMDalpha will be ignored with accelMDG on.\n" << endi;
	   if(opts.defined("accelMDTE"))
	       iout << iWARN << "accelMDTE will be ignored with accelMDG on.\n" << endi;
	   if(opts.defined("accelMDTalpha"))
	       iout << iWARN << "accelMDTalpha will be ignored with accelMDG on.\n" << endi;
       }
       else{
	   if(!opts.defined("accelMDE") || !opts.defined("accelMDalpha"))
	       NAMD_die("accelMDE and accelMDalpha are required for accelMD with accelMDG off");

	   if(accelMDdual && (!opts.defined("accelMDTE") || !opts.defined("accelMDTalpha"))){
	       NAMD_die("accelMDTE and accelMDTalpha are required for accelMDdual with accelMDG off");
	   }
       }
   }

   //  Set the default value for the maximum movement parameter
   //  for minimization
   if (minimizeOn && (maximumMove == 0.0))
   {
      maximumMove = 0.75 * pairlistDist/stepsPerCycle;
   }
   if (adaptTempOn) {
     if (!adaptTempRescale && !adaptTempLangevin)
        NAMD_die("Adaptive tempering needs to be coupled to either the Langevin thermostat or velocity rescaling.");
     if (opts.defined("adaptTempInFile") && (opts.defined("adaptTempTmin") ||
                                             opts.defined("adaptTempTmax") ||
                                             adaptTempBins != 0))
        NAMD_die("cannot simultaneously specify adaptTempInFile and any of {adaptTempTmin, adaptTempTmax,adaptTempBins} as these are read from the input file");
     if (!opts.defined("adaptTempInFile") && !(opts.defined("adaptTempTmin") &&
                                             opts.defined("adaptTempTmax") &&
                                             adaptTempBins != 0 ))
        NAMD_die("Need to specify either adaptTempInFile or all of {adaptTempTmin, adaptTempTmax,adaptTempBins} if adaptTempMD is on.");
   }

   langevinOnAtStartup = langevinOn;
   if (langevinOn) {
     if ( ! opts.defined("langevinDamping") ) langevinDamping = 0.0;
     if ( ! opts.defined("langevinHydrogen") ) langevinHydrogen = TRUE;
     if ( (opts.defined("langevinDamping") || opts.defined("langevinHydrogen"))
       && (opts.defined("langevinFile") || opts.defined("langevinCol")) )
       NAMD_die("To specify Langevin dynamics parameters, use either langevinDamping and langevinHydrogen or langevinFile and langevinCol.  Do not combine them.");
     if ( opts.defined("langevinHydrogen") && langevinDamping == 0.0 )
       NAMD_die("langevinHydrogen requires langevinDamping to be set.");
   }

   // BEGIN LA
   if (loweAndersenOn) {
       if (!opts.defined("loweAndersenRate")) loweAndersenRate = 100;
       if (!opts.defined("loweAndersenCutoff")) loweAndersenCutoff = 2.7;
   }
   // END LA

   // BKR - stochastic velocity rescaling
   if (stochRescaleOn) {
     if (langevinOn || loweAndersenOn || tCoupleOn ||
         opts.defined("rescaleFreq") || opts.defined("reassignFreq"))
       NAMD_die("Stochastic velocity rescaling is incompatible with other temperature control methods");
     // This is largely the same default used in GROMACS.
     if (!opts.defined("stochRescaleFreq")) stochRescaleFreq = stepsPerCycle;
   }

   if (opts.defined("rescaleFreq"))
   {
  if (!opts.defined("rescaleTemp"))
  {
    if (opts.defined("temperature"))
    {
      rescaleTemp = initialTemp;
    }
    else
    {
      NAMD_die("Must give a rescale temperature if rescaleFreq is defined");
    }
  }
   }
   else
   {
  rescaleFreq = -1;
  rescaleTemp = 0.0;
   }

   if (opts.defined("rescaleTemp"))
   {
  if (!opts.defined("rescaleFreq"))
  {
    NAMD_die("Must give a rescale freqency if rescaleTemp is given");
  }
   }

   if (opts.defined("reassignFreq"))
   {
  if (!opts.defined("reassignTemp"))
  {
    if (opts.defined("temperature"))
    {
      reassignTemp = initialTemp;
    }
    else
    {
      NAMD_die("Must give a reassign temperature if reassignFreq is defined");
    }
  }
   }
   else
   {
  reassignFreq = -1;
  reassignTemp = 0.0;
   }

   if (opts.defined("reassignTemp"))
   {
  if (!opts.defined("reassignFreq"))
  {
    NAMD_die("Must give a reassignment freqency if reassignTemp is given");
  }
   }

   if (opts.defined("reassignIncr"))
   {
  if (!opts.defined("reassignFreq"))
  {
    NAMD_die("Must give a reassignment freqency if reassignIncr is given");
  }
   }
   else
   {
  reassignIncr = 0.0;
   }

   if (opts.defined("reassignHold"))
   {
  if (!opts.defined("reassignIncr"))
  {
    NAMD_die("Must give a reassignment increment if reassignHold is given");
  }
   }
   else
   {
  reassignHold = 0.0;
   }

   if (!opts.defined("seed"))
   {
      randomSeed = (unsigned int) time(NULL) + 31530001 * CmiMyPartition();
   }

   // REST2
   if (opts.defined("soluteScaling")) {
     // Parameters soluteScalingFactorCharge and soluteScalingFactorVdw
     // allow independent scaling of electrostatics and van der Waals.
     // Initialize with soluteScalingFactor if either is not already set.
     if ( ! opts.defined("soluteScalingFactorCharge") ) {
       soluteScalingFactorCharge = soluteScalingFactor;
     }
     if ( ! opts.defined("soluteScalingFactorVdw") ) {
       soluteScalingFactorVdw = soluteScalingFactor;
     }
   }

//fepb
   alchFepOnAtStartup = alchFepOn = FALSE;
   alchThermIntOnAtStartup = alchThermIntOn = FALSE;
   alchOnAtStartup = alchOn;

   if (alchOn) {
     if (martiniSwitching) {
       iout << iWARN << "Martini switching disabled for alchemical "
         "interactions.\n" << endi;
     }

     if (!opts.defined("alchType")) {
       NAMD_die("Must define type of alchemical simulation: fep or ti\n");
     }
     else {
       opts.get("alchType",s);
       if (!strcasecmp(s, "fep")) {
         alchFepOnAtStartup = alchFepOn = TRUE;
       }
       else if (!strcasecmp(s, "ti")) {
         alchThermIntOnAtStartup = alchThermIntOn = TRUE;
       }
       else {
         NAMD_die("Unknown type of alchemical simulation; choices are fep or ti\n");
       }
     }

     if	     (rescaleFreq > 0) 	alchTemp = rescaleTemp;
     else if (reassignFreq > 0)	alchTemp = reassignTemp;
     else if (langevinOn) 	alchTemp = langevinTemp;
     else if (stochRescaleOn)	alchTemp = stochRescaleTemp;
     else if (tCoupleOn) 	alchTemp = tCoupleTemp;
     else NAMD_die("Alchemical FEP can be performed only in constant temperature simulations\n");

     if (reassignFreq > 0 && reassignIncr != 0)
	NAMD_die("reassignIncr cannot be used in alchemical simulations\n");

     if (alchLambda < 0.0 || alchLambda > 1.0)
       NAMD_die("alchLambda values should be in the range [0.0, 1.0]\n");

     if (alchVdwLambdaEnd > 1.0)
        NAMD_die("Gosh tiny Elvis, you kicked soft-core in the van der Waals! alchVdwLambdaEnd should be in the range [0.0, 1.0]\n");

     if (alchBondLambdaEnd > 1.0)
       NAMD_die("alchBondLambdaEnd should be in the range [0.0, 1.0]\n");

     if (alchElecLambdaStart > 1.0)
        NAMD_die("alchElecLambdaStart should be in the range [0.0, 1.0]\n");

     if (alchWCAOn) {
       if (alchRepLambdaEnd > 1.0)
         NAMD_die("alchRepLambdaEnd should be in the range [0.0, 1.0]\n");
       if (alchVdwLambdaEnd < alchRepLambdaEnd)
         NAMD_die("alchVdwLambdaEnd should be greater than alchRepLambdaEnd\n");
       if (alchVdwShiftCoeff > 0.0) {
         iout << iWARN << "alchVdwShiftCoeff is non-zero but not used when WCA"
              << " is active. Setting it to zero now.\n" << endi;
         alchVdwShiftCoeff = 0.0;
       }
       if (alchThermIntOn) {
         NAMD_die("alchWCA is not currently compatible with TI");
       }
#ifdef NAMD_CUDA
       NAMD_die("alchWCA is not currently available with CUDA");
#endif
     }

     if (alchFepOn) {
       if (alchLambda2 < 0.0 || alchLambda2 > 1.0)
         NAMD_die("alchLambda2 values should be in the range [0.0, 1.0]\n");

       setupIDWS(); // setup IDWS if it was activated.

       if (!opts.defined("alchoutfile")) {
         strcpy(alchOutFile, outputFilename);
         strcat(alchOutFile, ".fep");
       }

       if (!opts.defined("alchLambda") || !opts.defined("alchLambda2")) {
         NAMD_die("alchFepOn is on, but alchLambda or alchLambda2 is not set.");
       }
     }
     else if (alchThermIntOn) {
       // alchLambda2 is only needed for nonequilibrium switching
       if (alchLambdaFreq && (alchLambda2 < 0.0 || alchLambda2 > 1.0))
         NAMD_die("alchLambda2 values should be in the range [0.0, 1.0]\n");

       if (!opts.defined("alchoutfile")) {
         strcpy(alchOutFile, outputFilename);
         strcat(alchOutFile, ".ti");
       }
     }
   }

//fepe

   if ( alchOn && alchFepOn && alchThermIntOn )
     NAMD_die("Sorry, combined TI and FEP is not implemented.\n");
   if ( alchOn && lesOn )
     NAMD_die("Sorry, combined LES with FEP or TI is not implemented.\n");
   if ( alchOn && alchThermIntOn && lesOn )
     NAMD_die("Sorry, combined LES and TI is not implemented.\n");
   if ( alchWCAOn && !alchOn ) {
     iout << iWARN << "Alchemical WCA decomposition was requested but \
       alchemical free energy calculation is not active. Setting \
       alchWCA to off.\n" << endi;
     alchWCAOn = FALSE;
   }
   if ( alchDecouple && !alchOn ) {
         iout << iWARN << "Alchemical decoupling was requested but \
           alchemical free energy calculation is not active. Setting \
           alchDecouple to off.\n" << endi;
         alchDecouple = FALSE;
   }
   if ( alchBondDecouple && !alchOn ) {
         iout << iWARN << "Alchemical bond decoupling was requested but \
           alchemical free energy calculation is not active. Setting \
           alchBondDecouple to off.\n" << endi;
         alchBondDecouple = FALSE;
   }

   if ( lesOn && ( lesFactor < 1 || lesFactor > 255 ) ) {
     NAMD_die("lesFactor must be positive and less than 256");
   }
   if ((pairInteractionOn && alchOn) || (pairInteractionOn && lesOn))
     NAMD_die("Sorry, pair interactions may not be calculated when LES, FEP or TI is enabled.");

   // Drude model
   if (drudeOn) {
     if ( ! langevinOn ) {
       NAMD_die("Drude model requires use of Langevin thermostat.");
     }
     if ( ! opts.defined("drudeDamping")) {
       drudeDamping = langevinDamping;
       iout << iWARN << "Undefined 'drudeDamping' will be set to "
         "value of 'langevinDamping'\n" << endi;
     }
   }

   //  Set up load balancing variables
   if (opts.defined("ldBalancer")) {
     if (strcasecmp(loadBalancer, "none") == 0)
       ldBalancer = LDBAL_NONE;
     else if (strcasecmp(loadBalancer, "hybrid") == 0)
       ldBalancer = LDBAL_HYBRID;
     else
       NAMD_die("Unknown ldBalancer selected");
   } else {
     ldBalancer = LDBAL_CENTRALIZED;
#ifdef MEM_OPT_VERSION
     if ( CkNumPes() > 1400 ) ldBalancer = LDBAL_HYBRID;
#endif
   }

   if (opts.defined("ldbStrategy")) {
     //  Assign the load balancing strategy
     if (strcasecmp(loadStrategy, "comprehensive") == 0)
       ldbStrategy = LDBSTRAT_COMPREHENSIVE;
     else if (strcasecmp(loadStrategy, "refineonly") == 0)
       ldbStrategy = LDBSTRAT_REFINEONLY;
     else if (strcasecmp(loadStrategy, "old") == 0)
       ldbStrategy = LDBSTRAT_OLD;
     else
       NAMD_die("Unknown ldbStrategy selected");
   } else {
     ldbStrategy = LDBSTRAT_DEFAULT;
   }

  if (!opts.defined("ldbPeriod")) {
    ldbPeriod=200*stepsPerCycle;
  }

  //  Set default values
  if (!opts.defined("firstLdbStep")) {
    firstLdbStep=5*stepsPerCycle;
  }

  if (ldbPeriod <= firstLdbStep) {
    NAMD_die("ldbPeriod must greater than firstLdbStep.");
  }

  if (!opts.defined("lastLdbStep")) {
    lastLdbStep = -1;
  }

  if (!opts.defined("hybridGroupSize")) {
    hybridGroupSize = 512;
  }
  if ( hybridGroupSize < CkNumPes() ) {
    // match load balancer boundaries to physical nodes if possible
    int groupsize = hybridGroupSize;
    int *rpelist;
    int nodesize;
    CmiGetPesOnPhysicalNode(CmiPhysicalNodeID(0), &rpelist, &nodesize);
    if ( CkNumPes() % nodesize ) nodesize = CmiNodeSize(CmiNodeOf(0));
    if ( CkNumPes() % nodesize ) nodesize = 1;
    groupsize += nodesize - 1;
    while ( 2 * groupsize > CkNumPes() ) --groupsize;
    if ( groupsize < nodesize ) groupsize = nodesize;
    while ( groupsize % nodesize ) --groupsize;
    while ( groupsize && CkNumPes() % groupsize ) groupsize -= nodesize;
    if ( 2 * groupsize < hybridGroupSize ) {
      groupsize += nodesize;
      while ( CkNumPes() % groupsize ) groupsize += nodesize;
    }
    if ( 2 * groupsize <= CkNumPes() ) hybridGroupSize = groupsize;
  }

  // tracing will be done if trace is available and user says +traceOff
  // in that case we set nice values for some functions
  bool specialTracing = traceAvailable() && (traceIsOn() == 0);

  if(!opts.defined("traceStartStep")) {
    traceStartStep = 4 * firstLdbStep + 2 * ldbPeriod;
  }
  if(!opts.defined("numTraceSteps")) {
    numTraceSteps = 100;
  }

  if(specialTracing) {
    if (!opts.defined("firstLdbStep")) firstLdbStep = 20;
    if (!opts.defined("ldbPeriod")) ldbPeriod = 100;

    if(!opts.defined("traceStartStep")) {
      traceStartStep = 4 * firstLdbStep + 2 * ldbPeriod; // 380
    }

    if(!opts.defined("numTraceSteps")) {
      numTraceSteps = 80;
    }
  }

#ifdef MEASURE_NAMD_WITH_PAPI
  if(papiMeasure){
	  if(!opts.defined("papiMeasureStartStep")) {
		  papiMeasureStartStep = 3 * firstLdbStep;
	  }
	  if(!opts.defined("numPapiMeasureSteps")) {
		  numPapiMeasureSteps = 8; //including two pme steps
	  }
  }
#endif

  if(simulateInitialMapping) {
	  if(!opts.defined("simulatedPEs")){
		  simulatedPEs = CkNumPes();
	  }
	  if(!opts.defined("simulatedNodeSize")){
		  simulatedNodeSize = CkMyNodeSize();
	  }
  }

#ifdef MEM_OPT_VERSION
  //Some constraints on the values of load balancing parameters.
  //The reason is related to communication schemes used in sending proxy
  //data. If the step immediately after the load balancing is not a step
  //for atom migration, then it's possible there are some necessary information
  // missing inside the ProxyPatch which will crash the program. Therefore,
  // It's better that the step immediately after the load balancing be a step
  // for atom migration so that the some overhead in Proxy msgs are removed.
  // --Chao Mei
  if(ldbPeriod%stepsPerCycle!=0 || firstLdbStep%stepsPerCycle!=0) {
      iout << iWARN << "In memory optimized version, the ldbPeriod parameter or firstLdbStep parameter is better set to be a multiple of stepsPerCycle parameter!\n";
  }
#endif

   if (N < firstTimestep) { N = firstTimestep; }

   if ( (firstTimestep%stepsPerCycle) != 0)
   {
  NAMD_die("First timestep must be a multiple of stepsPerCycle!!");
   }

   //  Make sure only one full electrostatics algorithm is selected
   {
     int i = 0;
     if ( FMAOn ) ++i;
     if ( PMEOn ) ++i;
     if ( MSMOn ) ++i;
     if ( FMMOn ) ++i;
     if ( fullDirectOn ) ++i;
     if ( i > 1 )
	NAMD_die("More than one full electrostatics algorithm selected!!!");
   }

   if (!opts.defined("ldbBackgroundScaling")) {
     ldbBackgroundScaling = 1.0;
   }
   if (!opts.defined("ldbPMEBackgroundScaling")) {
     ldbPMEBackgroundScaling = ldbBackgroundScaling;
   }
   if (!opts.defined("ldbHomeBackgroundScaling")) {
     ldbHomeBackgroundScaling = ldbBackgroundScaling;
   }

   //  Check on PME parameters
   if (PMEOn) {  // idiot checking
     if ( lattice.volume() == 0. ) {
	NAMD_die("PME requires periodic boundary conditions.");
     }
     if ( PMEGridSpacing == 0. ) {
       if ( PMEGridSizeX * PMEGridSizeY * PMEGridSizeZ == 0 )
	 NAMD_die("Either PMEGridSpacing or PMEGridSizeX, PMEGridSizeY, and PMEGridSizeZ must be specified.");
       else PMEGridSpacing = 1.5;  // only exit in very bad cases
     }
#ifndef TEST_PME_GRID
     for ( int idim = 0; idim < 3; ++idim ) {
        int *gridSize;
        BigReal cellLength;
        const char *direction;
        switch ( idim ) {
        case 0:  direction = "X";
           gridSize = &PMEGridSizeX;  cellLength = lattice.a().length();
           break;
        case 1:  direction = "Y";
           gridSize = &PMEGridSizeY;  cellLength = lattice.b().length();
           break;
        case 2:  direction = "Z";
           gridSize = &PMEGridSizeZ;  cellLength = lattice.c().length();
           break;
        }
	int minSize = (int) ceil(cellLength/PMEGridSpacing);
#else
  for ( int minSize = 1; minSize < 300; ++minSize ) {
#endif
	int bestSize = 10 * (minSize + 10);  // make sure it's big
	int max2, max3, ts;
	for ( max2=2, ts=1; ts < minSize; ++max2 ) ts *= 2;
	for ( max3=2, ts=1; ts < minSize; ++max3 ) ts *= 3;
	int max5 = 2;
	int max7 = 1;
	int max11 = 1;
	for ( int i2 = 0; i2 <= max2; ++i2 ) {
	for ( int i3 = 0; i3 <= max3; ++i3 ) {
	for ( int i5 = 0; i5 <= max5; ++i5 ) {
	for ( int i7 = 0; i7 <= max7; ++i7 ) {
	for ( int i11 = 0; i11 <= max11; ++i11 ) {
	   if ( i5 + i7 + i11 > i2 ) continue;
           int testSize = 2;  // must be even
	   for ( int j2 = 0; j2 < i2; ++j2 ) testSize *= 2;
	   if ( testSize > bestSize ) continue;
	   for ( int j3 = 0; j3 < i3; ++j3 ) testSize *= 3;
	   if ( testSize > bestSize ) continue;
	   for ( int j5 = 0; j5 < i5; ++j5 ) testSize *= 5;
	   if ( testSize > bestSize ) continue;
	   for ( int j7 = 0; j7 < i7; ++j7 ) testSize *= 7;
	   if ( testSize > bestSize ) continue;
	   for ( int j11 = 0; j11 < i11; ++j11 ) testSize *= 11;
	   if ( testSize > bestSize ) continue;
	   if ( testSize >= minSize ) bestSize = testSize;
        } } } } }
#ifdef TEST_PME_GRID
  iout << minSize << " " << bestSize << "\n" << endi;
#else
	if ( ! *gridSize ) {   // set it
	   *gridSize = bestSize;
        }
	if ( *gridSize * PMEGridSpacing < cellLength ) {
	   char errmsg[512];
	   sprintf(errmsg, "PMEGridSize%s %d is too small for cell length %f and PMEGridSpacing %f\n",
		direction, *gridSize, cellLength, PMEGridSpacing);
	   NAMD_die(errmsg);
	}
#endif
     }
     if ( PMEGridSizeX < 5 ) {
	NAMD_die("PMEGridSizeX (number of grid points) is very small.");
     }
     if ( PMEGridSizeY < 5 ) {
	NAMD_die("PMEGridSizeY (number of grid points) is very small.");
     }
     if ( PMEGridSizeZ < 5 ) {
	NAMD_die("PMEGridSizeZ (number of grid points) is very small.");
     }
     BigReal tolerance = PMETolerance;
     BigReal ewaldcof = 1.0;
     while ( erfc(ewaldcof*cutoff)/cutoff >= tolerance ) ewaldcof *= 2.0;
     BigReal ewaldcof_lo = 0.;
     BigReal ewaldcof_hi = ewaldcof;
     for ( int i = 0; i < 100; ++i ) {
       ewaldcof = 0.5 * ( ewaldcof_lo + ewaldcof_hi );
       if ( erfc(ewaldcof*cutoff)/cutoff >= tolerance ) {
         ewaldcof_lo = ewaldcof;
       } else {
         ewaldcof_hi = ewaldcof;
       }
     }
     PMEEwaldCoefficient = ewaldcof;

#ifdef NAMD_CUDA
     bool one_device_per_node = deviceCUDA->one_device_per_node();  // only checks node 0
     if ( ! opts.defined("PMEOffload") ) {
       PMEOffload = ( (PMEInterpOrder > 4) && one_device_per_node );
       if ( PMEOffload ) iout << iINFO << "Enabling PMEOffload because PMEInterpOrder > 4.\n" << endi;
     } else if ( PMEOffload && ! one_device_per_node ) {
       PMEOffload = 0;
       iout << iWARN << "Disabling PMEOffload because multiple CUDA devices per process are not supported.\n" << endi;
     }
     if ( usePMECUDA && ! ( useCUDA2 || one_device_per_node ) ) {
       usePMECUDA = 0;
       iout << iWARN << "Disabling usePMECUDA because multiple CUDA devices per process requires useCUDA2.\n" << endi;
     }
#else
     PMEOffload = 0;
#endif
   } else {  // initialize anyway
     useDPME = 0;
     PMEGridSizeX = 0;
     PMEGridSizeY = 0;
     PMEGridSizeZ = 0;
     PMEGridSpacing = 1000.;
     PMEEwaldCoefficient = 0;
     PMEOffload = 0;
   }

   //  Take care of initializing FMA values to something if FMA is not
   //  active
   if (!FMAOn)
   {
     FMALevels = 0;
     FMAMp = 0;
     FMAFFTOn = FALSE;
     FMAFFTBlock = 0;
   }
   else
   {
  // idiot checking: frm bug reported by Tom Bishop.
  // DPMTA requires: (#terms in fma)/(fft blocking factor) = integer.
  if (FMAFFTBlock != 4)
    NAMD_die("FMAFFTBlock: Block length must be 4 for short FFT's");
  if (FMAMp % FMAFFTBlock != 0)
    NAMD_die("FMAMp: multipole term must be multiple of block length (FMAFFTBlock)");
    }

   if ( (nonbondedFrequency > stepsPerCycle) || ( (stepsPerCycle % nonbondedFrequency) != 0) )
   {
     NAMD_die("stepsPerCycle must be a multiple of nonbondedFreq");
   }

   if (!LCPOOn && !GBISOn && !GBISserOn && !FMAOn && !PMEOn && !MSMOn && !fullDirectOn && !FMMOn)
   {
     fullElectFrequency = 0;
   }
   else
   {
     if (!opts.defined("fullElectFrequency"))
     {
       if (opts.defined("fmaFrequency")) {
         iout << iWARN << "The parameter fmaFrequency has been renamed fullElectFrequency.\n" << endi;
         fullElectFrequency = fmaFrequency;
       } else {
         iout << iWARN << "The parameter fullElectFrequency now defaults to nonbondedFreq (" << nonbondedFrequency << ") rather than stepsPerCycle.\n" << endi;
         fullElectFrequency = nonbondedFrequency;
       }
     }
     else
     {
       if (opts.defined("fmaFrequency")) {
         iout << iWARN << "Ignoring redundant parameter fmaFrequency in favor of fullElectFrequency.\n" << endi;
       }
       if ( (fullElectFrequency > stepsPerCycle) || ( (stepsPerCycle % fullElectFrequency) != 0) )
       {
         NAMD_die("stepsPerCycle must be a multiple of fullElectFrequency");
       }
     }

     if ( (nonbondedFrequency > fullElectFrequency) || ( (fullElectFrequency % nonbondedFrequency) != 0) )
     {
       NAMD_die("fullElectFrequency must be a multiple of nonbondedFreq");
     }

     if (singleTopology && fullElectFrequency > 1) NAMD_die("Single topology free energy calculation discourages multiple timesteps to assure accuracy!");
     if (singleTopology && alchDecouple) NAMD_die("Single topology free energy calculation can NOT work with alchDecouple on");

      if (multigratorOn) {
        if ( (multigratorTemperatureFreq > multigratorPressureFreq) || ( (multigratorPressureFreq % multigratorTemperatureFreq) != 0) )
        {
          NAMD_die("multigratorTemperatureFreq must be a multiple of multigratorPressureFreq");
        }
        if ( (fullElectFrequency > multigratorTemperatureFreq) || ( (multigratorTemperatureFreq % fullElectFrequency) != 0) )
        {
          NAMD_die("fullElectFrequency must be a multiple of multigratorTemperatureFreq");
        }
        if (multigratorNoseHooverChainLength <= 2) {
          NAMD_die("multigratorNoseHooverChainLength must be greater than 2");
        }
      }

     if (!opts.defined("fmaTheta"))
     fmaTheta=0.715;  /* Suggested by Duke developers */
   }

   if ( lesOn && ( FMAOn || useDPME || fullDirectOn ) ) {
     NAMD_die("Sorry, LES is only implemented for PME full electrostatics.");
   }
   if ( alchFepOn && ( FMAOn || useDPME || fullDirectOn ) ) {
     NAMD_die("Sorry, FEP is only implemented for PME full electrostatics.");
   }
   if ( alchThermIntOn && ( FMAOn || useDPME || fullDirectOn ) ) {
     NAMD_die("Sorry, TI is only implemented for PME full electrostatics.");
   }
   if ( pairInteractionOn && FMAOn ) {
     NAMD_die("Sorry, pairInteraction not implemented for FMA.");
   }
   if ( pairInteractionOn && useDPME ) {
     NAMD_die("Sorry, pairInteraction not implemented for DPME.");
   }
   if ( pairInteractionOn && fullDirectOn ) {
     NAMD_die("Sorry, pairInteraction not implemented for full direct electrostatics.");
   }
   if ( ! pairInteractionOn ) {
     pairInteractionSelf = 0;
   }
   if ( pairInteractionOn && !pairInteractionSelf && !config->find("pairInteractionGroup2"))
     NAMD_die("pairInteractionGroup2 must be specified");

   if ( ! fixedAtomsOn ) {
     fixedAtomsForces = 0;
   }

   if ( gridforceOn || mgridforceOn ) {
     parse_mgrid_params(config);
   }

   if ( extraBondsOn ) {
     extraBondsCosAnglesSetByUser = ! ! config->find("extraBondsCosAngles");
   } else {
     extraBondsCosAnglesSetByUser = false;
   }

   if (!opts.defined("constraints"))
   {
     constraintExp = 0;
     constraintScaling = 1.0;

     //****** BEGIN selective restraints (X,Y,Z) changes
     selectConstraintsOn = FALSE;
     //****** END selective restraints (X,Y,Z) changes

     //****** BEGIN moving constraints changes
     movingConstraintsOn = FALSE;
     //****** END moving constraints changes
     //****** BEGIN rotating constraints changes
     rotConstraintsOn = FALSE;
    //****** END rotating constraints changes
   }
   //****** BEGIN rotating constraints changes
   else {
     if (rotConstraintsOn) {
       rotConsAxis = rotConsAxis.unit();
     }
   }
   if(opts.defined("rotConstraints")
      && opts.defined("movingConstraints")) {
     NAMD_die("Rotating and moving constraints are mutually exclusive!");
   }
   //****** END rotating constraints changes

   //****** BEGIN selective restraints (X,Y,Z) changes
   if(opts.defined("selectConstraints") && !opts.defined("selectConstrX")
      && !opts.defined("selectConstrY") && !opts.defined("selectConstrZ")) {
     NAMD_die("selectConstraints was specified, but no Cartesian components were defined!");
   }
   if (!opts.defined("selectConstraints")) {
       constrXOn = FALSE;
       constrYOn = FALSE;
       constrZOn = FALSE;
   }
   //****** END selective restraints (X,Y,Z) changes


   //****** BEGIN SMD constraints changes

   if (!opts.defined("SMD")) {
     SMDOn = FALSE;
   }

   if (SMDOn) {
     // normalize direction
     if (SMDDir.length2() == 0) {
       NAMD_die("SMD direction vector must be non-zero");
     }
     else {
       SMDDir = SMDDir.unit();
     }

     if (SMDOutputFreq > 0 && SMDOutputFreq < stepsPerCycle
	 || SMDOutputFreq % stepsPerCycle != 0) {
       NAMD_die("SMDOutputFreq must be a multiple of stepsPerCycle");
     }
   }

   //****** END SMD constraints changes

   if (!sphericalBCOn)
     {
  sphericalBCr1 = 0.0;
  sphericalBCk1 = 0.0;
  sphericalBCexp1 = 0;
  sphericalBCr2 = 0.0;
  sphericalBCk2 = 0.0;
  sphericalBCexp2 = 0;
   }
   else if (!opts.defined("sphericalBCr2"))
   {
      sphericalBCr2 = -1.0;
      sphericalBCk2 = 0.0;
      sphericalBCexp2 = 0;
   }

   if (!cylindricalBCOn)
   {
    cylindricalBCr1 = 0.0;
    cylindricalBCk1 = 0.0;
    cylindricalBCexp1 = 0;
    cylindricalBCr2 = 0.0;
    cylindricalBCk2 = 0.0;
    cylindricalBCexp2 = 0;
    cylindricalBCl1 = 0.0;
    cylindricalBCl2 = 0.0;
   }
   else if (!opts.defined("cylindricalBCr2"))
   {
    cylindricalBCr2 = -1.0;
    cylindricalBCk2 = 0.0;
    cylindricalBCexp2 = 0;
    cylindricalBCl2 = 0.0;
   }

   if (!eFieldOn)
   {
        eField.x = 0.0;
        eField.y = 0.0;
        eField.z = 0.0;
	eFieldFreq = 0.0;
	eFieldPhase = 0.0;
   }
   else
   {
   	if (!opts.defined("eFieldFreq")) eFieldFreq = 0.0;
	if (!opts.defined("eFieldPhase")) eFieldPhase = 0.0;
   }

   if (!stirOn)
   {
     stirFilename[0] = STRINGNULL;
     stirStartingTheta = 0.0;
     stirVel = 0.0;
     stirK = 0.0;
     stirAxis.x = 0.0;
     stirAxis.y = 0.0;
     stirAxis.z = 0.0;
     stirPivot.x = 0.0;
     stirPivot.y = 0.0;
     stirPivot.z = 0.0;
   }

   if (!opts.defined("langevin"))
   {
  langevinTemp = 0.0;
   }

   // BEGIN LA
   if (!opts.defined("loweAndersen"))
   {
       loweAndersenTemp = 0.0;
   }
   // END LA

   if (!opts.defined("tcouple"))
   {
  tCoupleTemp = 0.0;
   }

   if (HydrogenBonds)
   {
     if (daCutoffDist > pairlistDist)
       NAMD_die("Hydrogen bond cutoff distance must be <= pairlist distance");
   }

   // If we're doing pair interaction, set
   // outputEnergies to 1 to make NAMD not die (the other nonbonded code paths
   // aren't defined when these options are enabled), and set nonbondedFreq to
   // 1 to avoid getting erroneous output.  Warn the user of what we're doing.
   if (pairInteractionOn) {
	   if (outputEnergies != 1) {
		   iout << iWARN << "Setting outputEnergies to 1 due to\n";
		   iout << iWARN << "pairInteraction calculations\n" << endi;
		   outputEnergies  = 1;
	   }
   }
   if (pairInteractionOn || pressureProfileOn) {
	   if (nonbondedFrequency != 1) {
		   iout << iWARN << "Setting nonbondedFreq to 1 due to\n";
		   iout << iWARN << "pairInteraction or pressure profile calculations\n" << endi;
	   }
   }

   // print timing at a reasonable interval by default
   if (!opts.defined("outputTiming"))
   {
      outputTiming = firstLdbStep;
      int ot2 = 10 * outputEnergies;
      if ( outputTiming < ot2 ) outputTiming = ot2;
   }

    // Checks if a secondary process was added in the configuration, and sets
    // the appropriated variable
    if(qmForcesOn){

        if (opts.defined("QMSecProc")){
            qmSecProcOn = true;
        }
        else {
            qmSecProcOn = false;
        }

        if (opts.defined("qmPrepProc")){
            qmPrepProcOn = true;
        }
        else {
            qmPrepProcOn = false;
        }

        if (opts.defined("QMParamPDB")){
            qmParamPDBDefined = true;
        }
        else {
            qmParamPDBDefined = false;
        }

        if (opts.defined("QMBondColumn")){
            qmBondOn = true;
        }
        else {
            qmBondOn = false;
        }

        if ( strcasecmp(qmSoftware,"orca") != 0 &&
             strcasecmp(qmSoftware,"mopac") != 0 &&
             strcasecmp(qmSoftware,"custom") != 0 ) {
            NAMD_die("Available QM software options are \'mopac\', \'orca\', or \'custom\'.");
        }
        else {
            if ( strcasecmp(qmSoftware,"orca") == 0 )
                qmFormat = QMFormatORCA;
            if ( strcasecmp(qmSoftware,"mopac") == 0 )
                qmFormat = QMFormatMOPAC;
            if ( strcasecmp(qmSoftware,"custom") == 0 )
                qmFormat = QMFormatUSR;

            if (qmFormat == QMFormatORCA || qmFormat == QMFormatMOPAC) {

                if (! opts.defined("QMConfigLine"))
                    NAMD_die("If the selected QM software is \'mopac\' or \'orca\'\
, QMConfigLine needs to be defined.");

            }
        }

        qmChrgMode = QMCHRGMULLIKEN;
        if (opts.defined("QMChargeMode")) {
            if ( strcasecmp(qmChrgModeS,"none") != 0 &&
                 strcasecmp(qmChrgModeS,"mulliken") != 0 &&
                 strcasecmp(qmChrgModeS,"chelpg") != 0) {
                NAMD_die("Available charge options are \'none\', \'mulliken\' or \'chelpg\'.");
            }
            else {
                if ( strcasecmp(qmChrgModeS,"none") == 0 )
                    qmChrgMode = QMCHRGNONE;
                if ( strcasecmp(qmChrgModeS,"mulliken") == 0 )
                    qmChrgMode = QMCHRGMULLIKEN;
                if ( strcasecmp(qmChrgModeS,"chelpg") == 0 )
                    qmChrgMode = QMCHRGCHELPG;
            }
        }

        if (qmFormat == QMFormatMOPAC && qmChrgMode == QMCHRGCHELPG)
            NAMD_die("Available charge options for MOPAC are \'none\' and \'mulliken\'.");

        if (qmFormat == QMFormatUSR && qmChrgMode == QMCHRGCHELPG)
            NAMD_die("Available charge options for MOPAC are \'none\' and \'mulliken\'.");

        if (qmBondOn && (opts.defined("QMBondValueType"))) {
            if ( strcasecmp(qmBondValueTypeS,"len") != 0 &&
                strcasecmp(qmBondValueTypeS,"ratio") != 0 ) {
                NAMD_die("Available QM bond value type options are \'len\' or \'ratio\'.");
            }
            else {
    //         #define QMLENTYPE 1
    //         #define QMRATIOTYPE 2
                if ( strcasecmp(qmBondValueTypeS,"len") == 0 )
                    qmBondValType = 1;
                if ( strcasecmp(qmBondValueTypeS,"ratio") == 0 )
                    qmBondValType = 2;
            }
        }
        else if (qmBondOn && ! (opts.defined("QMBondValueType")))
            qmBondValType = 1;

        if ( strcmp(qmColumn,"beta") != 0 &&
             strcmp(qmColumn,"occ") != 0 ) {
            NAMD_die("Available column options are \'beta\' and \'occ\'.");
        }

        if (qmBondOn) {
            if ( strcmp(qmBondColumn,"beta") != 0 &&
                 strcmp(qmBondColumn,"occ") != 0 ) {
                NAMD_die("Available column options are \'beta\' and \'occ\'.");
            }

            if (strcmp(qmBondColumn,qmColumn) == 0)
                NAMD_die("QM column and bond-column must be different!");
        }

        qmBondScheme = 1;
        if (opts.defined("QMBondScheme")) {
            if ( strcasecmp(qmBondSchemeS,"CS") == 0 )
                qmBondScheme = QMSCHEMECS;
            if ( strcasecmp(qmBondSchemeS,"RCD") == 0 )
                qmBondScheme = QMSCHEMERCD;
            if ( strcasecmp(qmBondSchemeS,"Z1") == 0 )
                qmBondScheme = QMSCHEMEZ1;
            if ( strcasecmp(qmBondSchemeS,"Z2") == 0 )
                qmBondScheme = QMSCHEMEZ2;
            if ( strcasecmp(qmBondSchemeS,"Z3") == 0 )
                qmBondScheme = QMSCHEMEZ3;
        }

//         #define QMPCSCHEMENONE 1
//         #define QMPCSCHEMEROUND 2
//         #define QMPCSCHEMEZERO 3
        qmPCScheme = 1;
        if (opts.defined("QMPointChargeScheme") && qmPCSwitchOn) {
            if ( strcasecmp(qmPCSchemeS,"none") == 0 )
                qmPCScheme = 1;

            if ( strcasecmp(qmPCSchemeS,"round") == 0 )
                qmPCScheme = 2;
            if ( strcasecmp(qmPCSchemeS,"zero") == 0 )
                qmPCScheme = 3;

            if ( qmPCScheme > 1 && ! qmPCSwitchOn)
                NAMD_die("QM Charge Schemes \'round\' or \'zero\' can only be applied with QMswitching set to \'on\'!");
        }

        // Redundant option to deprecate "qmNoPC" option.
        if (qmElecEmbed)
            qmNoPC = FALSE;

//         #define QMLSSMODEDIST 1
//         #define QMLSSMODECOM 2
        if (qmLSSOn) {

            if (qmNoPC)
                NAMD_die("QM Live Solvent Selection cannot be done with QMNoPntChrg set to \'on\'!") ;

            if (rigidBonds != RIGID_NONE)
                NAMD_die("QM Live Solvent Selection cannot be done with fixed bonds!") ;

            if (qmLSSFreq % qmPCSelFreq != 0)
                NAMD_die("Frequency of QM solvent update must be a multiple of frequency of point charge selection.");

            if (qmLSSFreq % stepsPerCycle != 0)
                NAMD_die("Frequency of QM solvent update must be a multiple of steps per cycle.");

            if (opts.defined("QMLSSMode") ) {
                if ( strcasecmp(qmLSSModeS,"dist") != 0 &&
                     strcasecmp(qmLSSModeS,"COM") != 0 ) {
                    NAMD_die("Available LSS mode options are \'dist\' and \'COM\'.");
                }
                if ( strcasecmp(qmLSSModeS,"dist") == 0 )
                    qmLSSMode = 1;
                else if ( strcasecmp(qmLSSModeS,"COM") == 0 )
                    qmLSSMode = 2;
            }
            else
                qmLSSMode = 1;
        }

//         #define QMPCSCALESHIFT 1
//         #define QMPCSCALESWITCH 2
        if (qmPCSwitchOn) {

            if (opts.defined("QMSwitchingType") ) {
                if ( strcasecmp(qmPCSwitchTypeS,"shift") != 0 &&
                     strcasecmp(qmPCSwitchTypeS,"switch") != 0 ) {
                    NAMD_die("Available scaling options are \'shift\' and \'switch\'.");
                }
                if ( strcasecmp(qmPCSwitchTypeS,"shift") == 0 )
                    qmPCSwitchType = 1;
                else if ( strcasecmp(qmPCSwitchTypeS,"switch") == 0 )
                    qmPCSwitchType = 2;
            }
            else
                qmPCSwitchType = 1;
        }

        if (qmNoPC && qmPCSelFreq > 1) {
            iout << iWARN << "QMPCStride being IGNORED since QMNoPntChrg is set to \'on\'!\n" << endi;
            qmPCSelFreq = 1;
        }

        if (qmNoPC && qmPCSwitchOn)
            NAMD_die("QM PC switching can only be applied with QMNoPntChrg set to \'off\'!");

//         if (qmNoPC && qmBondOn)
//             NAMD_die("QM-MM bonds can only be applied with QMNoPntChrg set to \'off\'!");

        if (qmPCSelFreq <= 0)
            NAMD_die("QMPCFreq can only be a positive number! For static point charge selection, see QMCutomPC.");

        if (qmCustomPCSel && qmNoPC)
            NAMD_die("QM Custom PC Selection is incompatible with QMNoPntChrg!");

//         if (qmCustomPCSel && qmPCSwitchOn)
//             NAMD_die("QM Custom PC Selection is incompatible with QMSwitching!");

        if (qmCustomPCSel && qmPCSelFreq > 1)
            NAMD_die("QM Custom PC Selection is incompatible with QMPCStride > 1!");

        if (qmCSMD && (! opts.defined("QMCSMDFile") ))
            NAMD_die("QM Conditional SMD is ON, but no CSMD configuration file was profided!");
    }

#ifdef NAMD_CUDA
    // Disable various CUDA kernels if they do not fully support
    // or are otherwise incompatible with simulation options.
    if ( useCUDAdisable ) {
      if ( drudeOn && useCUDA2 && (bondedCUDA & 0x0001) ) {
        // disable CUDA kernels for spring bonds
        bondedCUDA &= ~0x0001;
        iout << iWARN << "Disabling CUDA kernel for bonds due to incompatibility with Drude oscillators.\n";
      }
      if ( accelMDOn && (accelMDdihe || accelMDdual) && useCUDA2 && (bondedCUDA & (0x0004 | 0x0020)) ) {
        // disable CUDA kernels for dihedrals and crossterms
        bondedCUDA &= ~(0x0004 | 0x0020);
        iout << iWARN << "Disabling CUDA kernels for dihedrals and crossterms due to incompatibility with accelerated MD options.\n";
      }
    }
#endif
} // check_config()


void SimParameters::print_config(ParseOptions &opts, ConfigList *config, char *&cwd) {

   StringList *current; //  Pointer to config option list

   //  Now that we have read everything, print it out so that
   //  the user knows what is going on
   iout << iINFO << "SIMULATION PARAMETERS:\n";
   iout << iINFO << "TIMESTEP               " << dt << "\n" << endi;
   iout << iINFO << "NUMBER OF STEPS        " << N << "\n";
   iout << iINFO << "STEPS PER CYCLE        " << stepsPerCycle << "\n";
   iout << endi;

   if ( lattice.a_p() || lattice.b_p() || lattice.c_p() ) {
     if ( lattice.a_p() )
       iout << iINFO << "PERIODIC CELL BASIS 1  " << lattice.a() << "\n";
     if ( lattice.b_p() )
       iout << iINFO << "PERIODIC CELL BASIS 2  " << lattice.b() << "\n";
     if ( lattice.c_p() )
       iout << iINFO << "PERIODIC CELL BASIS 3  " << lattice.c() << "\n";
     iout << iINFO << "PERIODIC CELL CENTER   " << lattice.origin() << "\n";
     if (wrapWater) {
       iout << iINFO << "WRAPPING WATERS AROUND PERIODIC BOUNDARIES ON OUTPUT.\n";
     }
     if (wrapAll) {
       iout << iINFO << "WRAPPING ALL CLUSTERS AROUND PERIODIC BOUNDARIES ON OUTPUT.\n";
     }
     if (wrapNearest) {
       iout << iINFO << "WRAPPING TO IMAGE NEAREST TO PERIODIC CELL CENTER.\n";
     }
     iout << endi;
   }

   if ( CkNumPes() > 512 ) ldbUnloadOne = TRUE;
   if ( ldbUnloadOne || CkNumPes() > 128 ) ldbUnloadZero = TRUE;

   if (ldBalancer == LDBAL_NONE) {
     iout << iINFO << "LOAD BALANCER  None\n" << endi;
   } else {
     if (ldBalancer == LDBAL_CENTRALIZED) {
       iout << iINFO << "LOAD BALANCER  Centralized\n" << endi;
     } else if (ldBalancer == LDBAL_HYBRID) {
       iout << iINFO << "LOAD BALANCER  Hybrid\n" << endi;
     }

     if (ldbStrategy == LDBSTRAT_DEFAULT) {
       iout << iINFO << "LOAD BALANCING STRATEGY  New Load Balancers -- DEFAULT\n";
     } else if (ldbStrategy == LDBSTRAT_REFINEONLY) {
       iout << iINFO << "LOAD BALANCING STRATEGY  Refinement Only\n";
     } else if (ldbStrategy == LDBSTRAT_COMPREHENSIVE) {
       iout << iINFO << "LOAD BALANCING STRATEGY  Comprehensive\n";
     } else if (ldbStrategy == LDBSTRAT_OLD) {
       iout << iINFO << "LOAD BALANCING STRATEGY  Old Load Balancers\n";
     }

     iout << iINFO << "LDB PERIOD             " << ldbPeriod << " steps\n";
     iout << iINFO << "FIRST LDB TIMESTEP     " << firstLdbStep << "\n";
     if (ldBalancer == LDBAL_HYBRID)
       iout << iINFO << "HYBRIDLB GROUP SIZE     " << hybridGroupSize << "\n";
     iout << iINFO << "LAST LDB TIMESTEP     " << lastLdbStep << "\n";
     if ( ldbRelativeGrainsize > 0. )
       iout << iINFO << "LDB RELATIVE GRAINSIZE " << ldbRelativeGrainsize << "\n";
     iout << iINFO << "LDB BACKGROUND SCALING " << ldbBackgroundScaling << "\n";
     iout << iINFO << "HOM BACKGROUND SCALING " << ldbHomeBackgroundScaling << "\n";
     if ( PMEOn ) {
       iout << iINFO << "PME BACKGROUND SCALING "
				<< ldbPMEBackgroundScaling << "\n";
     if ( ldbUnloadPME )
     iout << iINFO << "REMOVING LOAD FROM PME NODES" << "\n";
     }
     if ( ldbUnloadZero ) iout << iINFO << "REMOVING LOAD FROM NODE 0\n";
     if ( ldbUnloadOne ) iout << iINFO << "REMOVING LOAD FROM NODE 1\n";
     if ( ldbUnloadOutputPEs ) iout << iINFO << "REMOVING LOAD FROM OUTPUT PES\n";
     iout << endi;
   }

   if ( ldbUnloadOne || CkNumPes() > 256 ) noPatchesOnOne = TRUE;
   if ( ldbUnloadZero || noPatchesOnOne ||
          CkNumPes() > 64 || ( IMDon && CkNumPes() > 8 ) ) {
     noPatchesOnZero = TRUE;
   }
   if ( noPatchesOnZero ) iout << iINFO << "REMOVING PATCHES FROM PROCESSOR 0\n";
   if ( noPatchesOnOne ) iout << iINFO << "REMOVING PATCHES FROM PROCESSOR 1\n";
   iout << endi;

#if defined(NAMD_CUDA) || defined(NAMD_MIC)
    maxSelfPart = maxPairPart = 1;
#endif

   if (ldBalancer == LDBAL_HYBRID) {
     iout << iINFO << "MAX SELF PARTITIONS    " << maxSelfPart << "\n"
          << iINFO << "MAX PAIR PARTITIONS    " << maxPairPart << "\n"
          << iINFO << "SELF PARTITION ATOMS   " << numAtomsSelf << "\n"
          << iINFO << "SELF2 PARTITION ATOMS   " << numAtomsSelf2 << "\n"
          << iINFO << "PAIR PARTITION ATOMS   " << numAtomsPair << "\n"
          << iINFO << "PAIR2 PARTITION ATOMS  " << numAtomsPair2 << "\n";
   }
   iout << iINFO << "MIN ATOMS PER PATCH    " << minAtomsPerPatch << "\n"
        << endi;

   if (initialTemp < 0)
   {
  current = config->find("velocities");

  if (current == NULL)
  {
    current = config->find("binvelocities");
  }

  iout << iINFO << "VELOCITY FILE          " << current->data << "\n";
   }
   else
   {
  iout << iINFO << "INITIAL TEMPERATURE    "
     << initialTemp << "\n";
   }
   iout << endi;

   iout << iINFO << "CENTER OF MASS MOVING INITIALLY? ";

   if (comMove)
   {
     iout << "YES\n";
   }
   else
   {
     iout << "NO\n";
   }
   iout << endi;

   if ( zeroMomentum ) {
     iout << iINFO << "REMOVING CENTER OF MASS DRIFT DURING SIMULATION";
     if ( zeroMomentumAlt ) iout << " (ALT METHOD)";
     iout << "\n" << endi;
   }

   iout << iINFO << "DIELECTRIC             "
      << dielectric << "\n";

   if ( nonbondedScaling != 1.0 )
   {
     iout << iINFO << "NONBONDED SCALING    " << nonbondedScaling << "\n" << endi;
   }
   iout << iINFO << "EXCLUDE                ";

   switch (exclude)
   {
     case NONE:
       iout << "NONE\n";
       break;
     case ONETWO:
       iout << "ONETWO\n";
       break;
     case ONETHREE:
       iout << "ONETHREE\n";
       break;
     case ONEFOUR:
       iout << "ONE-FOUR\n";
       break;
     default:
       iout << "SCALED ONE-FOUR\n";
       break;
   }
   iout << endi;

   if (exclude == SCALED14)
   {
     iout << iINFO << "1-4 ELECTROSTATICS SCALED BY " << scale14 << "\n";
     iout << iINFO << "MODIFIED 1-4 VDW PARAMETERS WILL BE USED\n" << endi;
   } else {
     iout << iWARN << "MODIFIED 1-4 VDW PARAMETERS WILL BE IGNORED\n" << endi;
   }

#ifdef SPEC_DISABLED_VERSION
   if (dcdFrequency > 0) {
     dcdFrequency = 0;
     iout << iWARN << "DCD TRAJECTORY OUTPUT IS DISABLED IN SPEC RELEASE\n";
   }
#endif

   if (dcdFrequency > 0)
   {
     iout << iINFO << "DCD FILENAME           "
        << dcdFilename << "\n";
     iout << iINFO << "DCD FREQUENCY          "
        << dcdFrequency << "\n";
     iout << iINFO << "DCD FIRST STEP         "
        << ( ((firstTimestep + dcdFrequency)/dcdFrequency)*dcdFrequency ) << "\n";
     if ( dcdUnitCell ) {
       iout << iINFO << "DCD FILE WILL CONTAIN UNIT CELL DATA\n";
     }
   }
   else
   {
     iout << iINFO << "NO DCD TRAJECTORY OUTPUT\n";
   }
   iout << endi;

   if (xstFrequency > 0)
   {
     iout << iINFO << "XST FILENAME           "
        << xstFilename << "\n";
     iout << iINFO << "XST FREQUENCY          "
        << xstFrequency << "\n";
   }
   else
   {
     iout << iINFO << "NO EXTENDED SYSTEM TRAJECTORY OUTPUT\n";
   }
   iout << endi;

   if (velDcdFrequency > 0)
   {
     iout << iINFO << "VELOCITY DCD FILENAME    "
        << velDcdFilename << "\n";
     iout << iINFO << "VELOCITY DCD FREQUENCY   "
        << velDcdFrequency << "\n";
     iout << iINFO << "VELOCITY DCD FIRST STEP  "
        << ( ((firstTimestep + velDcdFrequency)/velDcdFrequency)*velDcdFrequency ) << "\n";
   }
   else
   {
     iout << iINFO << "NO VELOCITY DCD OUTPUT\n";
   }
   iout << endi;

   if (forceDcdFrequency > 0)
   {
     iout << iINFO << "FORCE DCD FILENAME     "
        << forceDcdFilename << "\n";
     iout << iINFO << "FORCE DCD FREQUENCY    "
        << forceDcdFrequency << "\n";
     iout << iINFO << "FORCE DCD FIRST STEP   "
        << ( ((firstTimestep + forceDcdFrequency)/forceDcdFrequency)*forceDcdFrequency ) << "\n";
   }
   else
   {
     iout << iINFO << "NO FORCE DCD OUTPUT\n";
   }
   iout << endi;

   iout << iINFO << "OUTPUT FILENAME        "
      << outputFilename << "\n" << endi;
   if (binaryOutput)
   {
     iout << iINFO << "BINARY OUTPUT FILES WILL BE USED\n" << endi;
   }
#ifdef MEM_OPT_VERSION
    if(!binaryOutput){
	iout << iWARN <<"SINCE MEMORY OPTIMIZED VERSION IS USED, OUTPUT IN TEXT FORMAT IS DISABLED!\n" << endi;
	binaryOutput = TRUE;
    }
#endif

   if (! restartFrequency)
   {
     iout << iINFO << "NO RESTART FILE\n";
   }
   else
   {
     iout << iINFO << "RESTART FILENAME       "
        << restartFilename << "\n";
     iout << iINFO << "RESTART FREQUENCY      "
        << restartFrequency << "\n";
  if (restartSave) {
    iout << iINFO << "RESTART FILES WILL NOT BE OVERWRITTEN\n";
  }
  if (restartSaveDcd) {
    iout << iINFO << "DCD FILE WILL BE SPLIT WHEN RESTART FILES ARE WRITTEN\n";
  }

  if (binaryRestart)
  {
    iout << iINFO << "BINARY RESTART FILES WILL BE USED\n";
  }
   }
   iout << endi;

   if (switchingActive)
   {
      iout << iINFO << "SWITCHING ACTIVE\n";
      if ( vdwForceSwitching ) {
        iout << iINFO << "VDW FORCE SWITCHING ACTIVE\n";
      }
      if ( martiniSwitching ) {
        iout << iINFO << "MARTINI RESIDUE-BASED COARSE-GRAIN SWITCHING ACTIVE\n";
      }
      iout << iINFO << "SWITCHING ON           "
               << switchingDist << "\n";
      iout << iINFO << "SWITCHING OFF          "
               << cutoff << "\n";
   }
   else
   {
      iout << iINFO << "CUTOFF                 "
         << cutoff << "\n";
   }

   iout << iINFO << "PAIRLIST DISTANCE      " << pairlistDist << "\n";
   iout << iINFO << "PAIRLIST SHRINK RATE   " << pairlistShrink << "\n";
   iout << iINFO << "PAIRLIST GROW RATE     " << pairlistGrow << "\n";
   iout << iINFO << "PAIRLIST TRIGGER       " << pairlistTrigger << "\n";
   iout << iINFO << "PAIRLISTS PER CYCLE    " << pairlistsPerCycle << "\n";
   if ( outputPairlists )
     iout << iINFO << "PAIRLIST OUTPUT STEPS  " << outputPairlists << "\n";
   iout << endi;

   if ( pairlistMinProcs > 1 )
     iout << iINFO << "REQUIRING " << pairlistMinProcs << " PROCESSORS FOR PAIRLISTS\n";
   usePairlists = ( CkNumPes() >= pairlistMinProcs );

#ifdef OPENATOM_VERSION
if ( openatomOn )
{
  iout << iINFO << "OPENATOM QM/MM CAR-PARINELLO ACTIVE\n";
  iout << iINFO << "OPENATOM CONFIG FILE:  " << openatomConfig << "\n";
  iout << iINFO << "OPENATOM STRUCT FILE:  " << openatomStruct << "\n";
  iout << iINFO << "OPENATOM PDB FILE:     " << openatomPDB << "\n";
}
#endif // OPENATOM_VERSION

   // FB - FEP and TI are now dependent on pairlists - disallow usePairlists=0
   if ( (alchOn) && (!usePairlists)) {
     NAMD_die("Sorry, Alchemical simulations require pairlists to be enabled\n");
   }
#ifdef NAMD_CUDA
   if ( ! usePairlists ) {
     usePairlists = 1;
     iout << iINFO << "CUDA ACCELERATION REQUIRES PAIRLISTS\n";
   }
#endif

   iout << iINFO << "PAIRLISTS " << ( usePairlists ? "ENABLED" : "DISABLED" )
							<< "\n" << endi;

   iout << iINFO << "MARGIN                 " << margin << "\n";
   if ( margin > 4.0 ) {
      iout << iWARN << "MARGIN IS UNUSUALLY LARGE AND WILL LOWER PERFORMANCE\n";
      BigReal f = patchDimension/(patchDimension-margin);
      f *= f*f;
      iout << iWARN << "MARGIN INCREASED PATCH VOLUME BY A FACTOR OF " << f << "\n";
   }

   if ( splitPatch == SPLIT_PATCH_HYDROGEN ) {
      iout << iINFO << "HYDROGEN GROUP CUTOFF  " << hgroupCutoff << "\n";
   }

   iout << iINFO << "PATCH DIMENSION        "
            << patchDimension << "\n";

   iout << endi;

   if (outputEnergies != 1)
   {
      iout << iINFO << "ENERGY OUTPUT STEPS    "
         << outputEnergies << "\n";
      iout << endi;
   }

   if (mergeCrossterms) {
      iout << iINFO << "CROSSTERM ENERGY INCLUDED IN DIHEDRAL\n" << endi;
   }

   if (outputMomenta != 0)
   {
      iout << iINFO << "MOMENTUM OUTPUT STEPS  "
         << outputMomenta << "\n";
      iout << endi;
   }

   if (outputTiming != 0)
   {
      iout << iINFO << "TIMING OUTPUT STEPS    "
         << outputTiming << "\n";
      iout << endi;
   }

   if (outputCudaTiming != 0)
   {
      iout << iINFO << "CUDA TIMING OUTPUT STEPS    "
         << outputCudaTiming << "\n";
      iout << endi;
   }

   if (outputPressure != 0)
   {
      iout << iINFO << "PRESSURE OUTPUT STEPS  "
         << outputPressure << "\n";
      iout << endi;
   }

   if (fixedAtomsOn)
   {
      iout << iINFO << "FIXED ATOMS ACTIVE\n";
      if ( fixedAtomsForces )
	iout << iINFO << "FORCES BETWEEN FIXED ATOMS ARE CALCULATED\n";
      iout << endi;
   }

   if (constraintsOn)
   {
      iout << iINFO << "HARMONIC CONSTRAINTS ACTIVE\n";

      iout << iINFO << "HARMONIC CONS EXP      "
         << constraintExp << "\n";

      if (constraintScaling != 1.0) {
        iout << iINFO << "HARMONIC CONS SCALING  "
         << constraintScaling << "\n";
      }

      //****** BEGIN selective restraints (X,Y,Z) changes

      if (selectConstraintsOn) {
	iout << iINFO << "SELECTED CARTESIAN COMPONENTS OF HARMONIC RESTRAINTS ACTIVE\n";

        if (constrXOn)
	iout << iINFO << "RESTRAINING X-COMPONENTS OF CARTESIAN COORDINATES!\n";

        if (constrYOn)
	iout << iINFO << "RESTRAINING Y-COMPONENTS OF CARTESIAN COORDINATES!\n";

        if (constrZOn)
	iout << iINFO << "RESTRAINING Z-COMPONENTS OF CARTESIAN COORDINATES!\n";
      }
      //****** END selective restraints (X,Y,Z) changes

      if (sphericalConstraintsOn) {
	iout << iINFO << "SPHERICAL HARMONIC CONSTRAINTS ACTIVE\n";
	iout << iINFO << "RESTRAINING DISTANCE TO " << sphericalConstrCenter <<"\n";
      }
      iout << endi;

      //****** BEGIN moving constraints changes

      if (movingConstraintsOn) {
	iout << iINFO << "MOVING HARMONIC CONSTRAINTS ACTIVE\n";

	iout << iINFO << "MOVING CONSTRAINT VELOCITY    "
	     << movingConsVel << " ANGSTROM/TIMESTEP\n";

	iout << iINFO << "ALL CONSTRAINED ATOMS WILL MOVE\n";
      }
      //****** END moving constraints changes
      iout << endi;

      //****** BEGIN rotating constraints changes

      if (rotConstraintsOn) {
	iout << iINFO << "ROTATING HARMONIC CONSTRAINTS ACTIVE\n";

	iout << iINFO << "AXIS OF ROTATION    "
	     << rotConsAxis << "\n";

	iout << iINFO << "PIVOT OF ROTATION   "
	     << rotConsPivot << "\n";

	iout << iINFO << "ROTATING CONSTRAINT VELOCITY    "
	     << rotConsVel << " DEGREES/TIMESTEP\n";
      }
      iout << endi;
      //****** END rotating constraints changes
   }

   // moving drag
   if (movDragOn) {
     iout << iINFO << "MOVING DRAG ACTIVE.\n";

     iout << iINFO << "MOVING DRAG MAIN PDB FILE "
	  << movDragFile << "\n";

     iout << iINFO << "MOVING DRAG GLOBAL VELOCITY (A/step) "
	  << movDragGlobVel << "\n";

     iout << iINFO << "MOVING DRAG LINEAR VELOCITY FILE "
	  << movDragVelFile << "\n";

     iout << endi;
   }

   // rotating drag
   if (rotDragOn) {
     iout << iINFO << "ROTATING DRAG ACTIVE.\n";

     iout << iINFO << "ROTATING DRAG MAIN PDB FILE "
	  << rotDragFile << "\n";

     iout << iINFO << "ROTATING DRAG AXIS FILE "
	  << rotDragAxisFile << "\n";

     iout << iINFO << "ROTATING DRAG PIVOT POINT FILE "
	  << rotDragPivotFile << "\n";

     iout << iINFO << "ROTATING DRAG GLOBAL ANGULAR VELOCITY (deg/step) "
	  << rotDragGlobVel << "\n";

     iout << iINFO << "ROTATING DRAG ANGULAR VELOCITY FILE "
	  << rotDragVelFile << "\n";

     iout << endi;
   }


   // "constant" torque
   if (consTorqueOn) {
     iout << iINFO << "\"CONSTANT\" TORQUE ACTIVE.\n";

     iout << iINFO << "\"CONSTANT\" TORQUE MAIN PDB FILE "
	  << consTorqueFile << "\n";

     iout << iINFO << "\"CONSTANT\" TORQUE AXIS FILE "
	  << consTorqueAxisFile << "\n";

     iout << iINFO << "\"CONSTANT\" TORQUE PIVOT POINT FILE "
	  << consTorquePivotFile << "\n";

     iout << iINFO << "\"CONSTANT\" TORQUE GLOBAL VALUE (Kcal/(mol*A^2)) "
	  << consTorqueGlobVal << "\n";

     iout << iINFO << "\"CONSTANT\" TORQUE DACTORS FILE "
	  << consTorqueValFile << "\n";

     iout << endi;
   }

   if (mgridforceOn) {
     iout << iINFO << "GRID FORCE ACTIVE\n";
     iout << iINFO << " Please include this reference in published work using\n";
     iout << iINFO << " the Gridforce module of NAMD: David Wells, Volha Abramkina,\n";
     iout << iINFO << " and Aleksei Aksimentiev, J. Chem. Phys. 127:125101-10 (2007).\n";
     print_mgrid_params();
   }

   //****** BEGIN SMD constraints changes

   if (SMDOn) {
     iout << iINFO << "SMD ACTIVE\n";

     iout << iINFO << "SMD VELOCITY    "
	  << SMDVel << " ANGSTROM/TIMESTEP\n";

     iout << iINFO << "SMD DIRECTION   "
	  << SMDDir << "\n";

     iout << iINFO << "SMD K   "
          << SMDk << "\n";

     iout << iINFO << "SMD K2  "
          << SMDk2 << "\n";

     iout << iINFO << "SMD OUTPUT FREQUENCY   "
	  << SMDOutputFreq << " TIMESTEPS\n";

     iout << iINFO << "SMD FILE " << SMDFile << "\n";

     iout << endi;
   }

   //****** END SMD constraints changes

   if (TMDOn) {
     iout << iINFO << "TMD ACTIVE BETWEEN STEPS " << TMDFirstStep
          << " and " << TMDLastStep << "\n";
     iout << iINFO << "TMD K  " << TMDk << "\n";
     iout << iINFO << "TMD FILE  " << TMDFile << "\n";
     iout << iINFO << "TMD OUTPUT FREQUENCY  " << TMDOutputFreq << "\n";
     if (TMDInitialRMSD) {
       iout << iINFO << "TMD TARGET RMSD AT FIRST STEP  " << TMDInitialRMSD << "\n";
     } else {
       iout << iINFO << "TMD TARGET RMSD AT FIRST STEP COMPUTED FROM INITIAL COORDINATES\n";
     }
     iout << iINFO << "TMD TARGET RMSD AT FINAL STEP  " << TMDFinalRMSD << "\n";
     iout << endi;
   }

   if (symmetryOn) {
     if (symmetryLastStep == -1){
       iout << iINFO << "SYMMETRY RESTRAINTS ACTIVE BETWEEN STEPS " << symmetryFirstStep << " and " << "INFINITY" << "\n";
     }
     else{
       iout << iINFO << "SYMMETRY RESTRAINTS ACTIVE BETWEEN STEPS " << symmetryFirstStep << " and " << symmetryLastStep << "\n";
     }
    // iout << iINFO << "SYMMETRY FILE " << symmetryFile << "\n";

     current = config->find("symmetryFile");
     for ( ; current; current = current->next ) {
       iout << iINFO << "SYMMETRY FILE  " << current->data << "\n";
     }

     current = config->find("symmetryMatrixFile");
     for ( ; current; current = current->next ) {
      iout << iINFO << "SYMMETRY MATRIX FILE " << current->data << "\n";
     }
     iout << iINFO << "SYMMETRY FORCE CONSTANT " << symmetryk << "\n";
     if (symmetryScaleForces){
      iout << iINFO << "SYMMETRY SCALE FORCES ON\n";
     }
     iout << iINFO << "SYMMETRY FIRST FULL STEP " << symmetryFirstFullStep << "\n";
     if (symmetryLastFullStep == -1){
      iout << iINFO << "SYMMETRY LAST FULL STEP " << "INFINITY" << "\n";
      //iout << iINFO << "FULL SYMMETRY FORCE BETWEEN STEPS " << symmetryFirstFullStep << " and " << "INFINITY" << "\n";
     }
     else {
      iout << iINFO << "SYMMETRY LAST FULL STEP " << symmetryLastFullStep << "\n";
     // iout << iINFO << "FULL SYMMETRY FORCE BETWEEN STEPS " << symmetryFirstFullStep << " and " << symmetryLastFullStep << "\n";
     }

     iout << endi;
   }
//Modifications for alchemical fep
//  Alchemical FEP status

//   current = config->find("alchOutFile");
   if (alchFepOn)
   {
     iout << iINFO << "ALCHEMICAL FEP ON\n";
     iout << iINFO << "FEP CURRENT LAMBDA VALUE     "
          << alchLambda << "\n";
     iout << iINFO << "FEP COMPARISON LAMBDA VALUE  "
          << alchLambda2 << "\n";
     if (alchLambdaIDWS >= 0.) {
        iout << iINFO << "FEP ALTERNATE COMPARISON LAMBDA VALUE  "
          << alchLambdaIDWS << "\n";
     }
     if (alchLambdaFreq > 0) {
       iout << iINFO << "FEP CURRENT LAMBDA VALUE SET TO INCREASE IN EVERY  "
            << alchLambdaFreq << " STEPS\n";
     }
     if (!alchDecouple) {
       iout << iINFO << "FEP INTRA-ALCHEMICAL NON-BONDED INTERACTIONS WILL BE "
            << "DECOUPLED\n";
     }else{
       iout << iINFO << "FEP INTRA-ALCHEMICAL NON-BONDED INTERACTIONS WILL BE "
            << "RETAINED\n";
     }
     if (alchBondDecouple) {
       iout << iINFO << "FEP INTRA-ALCHEMICAL BONDED INTERACTIONS WILL BE "
            << "DECOUPLED\n";
     }else{
       iout << iINFO << "FEP INTRA-ALCHEMICAL BONDED INTERACTIONS WILL BE "
            << "RETAINED\n";
     }
     if (alchWCAOn) {
       iout << iINFO << "FEP WEEKS-CHANDLER-ANDERSEN (WCA) VDW DECOUPLING "
            << "ACTIVE\n";
     } else {
     iout << iINFO << "FEP VDW SHIFTING COEFFICIENT "
          << alchVdwShiftCoeff << "\n";
     }
     iout << iINFO << "FEP ELEC. ACTIVE FOR ANNIHILATED "
          << "PARTICLES BETWEEN LAMBDA = 0 AND LAMBDA = "
          << (1 - alchElecLambdaStart) << "\n";
     iout << iINFO << "FEP ELEC. ACTIVE FOR EXNIHILATED "
          << "PARTICLES BETWEEN LAMBDA = "
          << alchElecLambdaStart << " AND LAMBDA = 1\n";
     if (alchWCAOn) {
       iout << iINFO << "FEP VDW-REPU. ACTIVE FOR ANNIHILATED PARTICLES "
            << "BETWEEN LAMBDA = " << (1 - alchRepLambdaEnd) << " AND LAMBDA "
            << "= 1\n";
       iout << iINFO << "FEP VDW-REPU. ACTIVE FOR EXNIHILATED PARTICLES "
            << "BETWEEN LAMBDA = 0 AND LAMBDA " << alchRepLambdaEnd << "\n";
       iout << iINFO << "FEP VDW-ATTR. ACTIVE FOR ANNIHILATED PARTICLES "
            << "BETWEEN LAMBDA = " << (1 - alchVdwLambdaEnd) << " AND LAMBDA = "
            << (1 - alchRepLambdaEnd) << "\n";
       iout << iINFO << "FEP VDW-ATTR. ACTIVE FOR EXNIHILATED PARTICLES "
            << "BETWEEN LAMBDA = " << alchRepLambdaEnd << " AND LAMBDA = "
            << alchVdwLambdaEnd << "\n";
     } else {
     iout << iINFO << "FEP VDW ACTIVE FOR ANNIHILATED "
          << "PARTICLES BETWEEN LAMBDA = "
          << (1 - alchVdwLambdaEnd) << " AND LAMBDA = 1\n";
     iout << iINFO << "FEP VDW ACTIVE FOR EXNIHILATED "
          << "PARTICLES BETWEEN LAMBDA = 0 AND LAMBDA = "
          << alchVdwLambdaEnd << "\n";
     }
     iout << iINFO << "FEP BOND ACTIVE FOR ANNIHILATED "
          << "PARTICLES BETWEEN LAMBDA = "
          << (1 - alchBondLambdaEnd) << " AND LAMBDA = 1\n";
     iout << iINFO << "FEP BOND ACTIVE FOR EXNIHILATED "
          << "PARTICLES BETWEEN LAMBDA = 0 AND LAMBDA = "
          << alchBondLambdaEnd << "\n";
   }
//fepe

   if (alchThermIntOn)
   {
     iout << iINFO << "THERMODYNAMIC INTEGRATION (TI) ON\n";
     iout << iINFO << "TI LAMBDA VALUE     "
          << alchLambda << "\n";
     if (alchLambdaFreq > 0) {
       iout << iINFO << "TI COMPARISON LAMBDA VALUE  "
            << alchLambda2 << "\n";
       iout << iINFO << "TI CURRENT LAMBDA VALUE SET TO INCREASE IN EVERY  "
            << alchLambdaFreq << " STEPS\n";
     }
     if (!alchDecouple) {
       iout << iINFO << "TI INTRA-ALCHEMICAL NON-BONDED INTERACTIONS WILL BE "
            << "DECOUPLED\n";
     }else{
       iout << iINFO << "TI INTRA-ALCHEMICAL NON-BONDED INTERACTIONS WILL BE "
            << "RETAINED\n";
     }
     if (alchBondDecouple) {
       iout << iINFO << "TI INTRA-ALCHEMICAL BONDED INTERACTIONS WILL BE "
            << "DECOUPLED\n";
     }else{
       iout << iINFO << "TI INTRA-ALCHEMICAL BONDED INTERACTIONS WILL BE "
            << "RETAINED\n";
     }
     iout << iINFO << "TI VDW SHIFTING COEFFICIENT "
          << alchVdwShiftCoeff << "\n";
     iout << iINFO << "TI ELEC. ACTIVE FOR ANNIHILATED "
          << "PARTICLES BETWEEN LAMBDA = 0 AND LAMBDA = "
          << (1 - alchElecLambdaStart) << "\n";
     iout << iINFO << "TI ELEC. ACTIVE FOR EXNIHILATED "
          << "PARTICLES BETWEEN LAMBDA = "
          << alchElecLambdaStart << " AND LAMBDA = 1\n";
     iout << iINFO << "TI VDW ACTIVE FOR ANNIHILATED "
          << "PARTICLES BETWEEN LAMBDA = "
          << (1 - alchVdwLambdaEnd) << " AND LAMBDA = 1\n";
     iout << iINFO << "TI VDW ACTIVE FOR EXNIHILATED "
          << "PARTICLES BETWEEN LAMBDA = 0 AND LAMBDA = "
          << alchVdwLambdaEnd << "\n";
     iout << iINFO << "TI BOND ACTIVE FOR ANNIHILATED "
          << "PARTICLES BETWEEN LAMBDA = "
          << (1 - alchBondLambdaEnd) << " AND LAMBDA = 1\n";
     iout << iINFO << "TI BOND ACTIVE FOR EXNIHILATED "
          << "PARTICLES BETWEEN LAMBDA = 0 AND LAMBDA = "
          << alchBondLambdaEnd << "\n";
   }


   if ( lesOn ) {
     iout << iINFO << "LOCALLY ENHANCED SAMPLING ACTIVE\n";
     iout << iINFO << "LOCAL ENHANCEMENT FACTOR IS "
          << lesFactor << "\n";
     if ( lesReduceTemp ) iout << iINFO
       << "SCALING ENHANCED ATOM TEMPERATURE BY 1/" << lesFactor << "\n";
     if ( lesReduceMass ) iout << iINFO
       << "SCALING ENHANCED ATOM MASS BY 1/" << lesFactor << "\n";
   }

   if ( singleTopology ) {
     iout << iINFO << "SINGLE TOPOLOGY IS ON FOR RELATIVE FREE ENERGY CALCULATION\n";
       }

   // REST2
   if ( soluteScalingOn ) {
     iout << iINFO << "SOLUTE SCALING IS ACTIVE\n";
     if (soluteScalingFactorCharge != soluteScalingFactorVdw) {
       iout << iINFO << "SCALING FOR ELECTROSTATIC INTERACTIONS IS "
            << soluteScalingFactorCharge << "\n";
       iout << iINFO << "SCALING FOR VAN DER WAALS INTERACTIONS IS "
            << soluteScalingFactorVdw << "\n";
       iout << iINFO << "SCALING FOR BONDED INTERACTIONS IS "
            << soluteScalingFactor << "\n";
     }
     else {
       iout << iINFO << "SOLUTE SCALING FACTOR IS "
            << soluteScalingFactor << "\n";
     }
     if ( ! soluteScalingAll ) {
       iout << iINFO << "SOLUTE SCALING DISABLED FOR BONDS AND ANGLES\n";
     }
   }

   if ( pairInteractionOn ) {
     iout << iINFO << "PAIR INTERACTION CALCULATIONS ACTIVE\n";
     iout << iINFO << "USING FLAG " << pairInteractionGroup1
          << " FOR GROUP 1\n";
     if (pairInteractionSelf) {
       iout << iINFO << "COMPUTING ONLY SELF INTERACTIONS FOR GROUP 1 ATOMS\n";
     } else {
       iout << iINFO << "USING FLAG " << pairInteractionGroup2
            << " FOR GROUP 2\n";
     }
   }

   if (consForceOn) {
     iout << iINFO << "CONSTANT FORCE ACTIVE\n";
     if ( consForceScaling != 1.0 ) {
       iout << iINFO << "CONSTANT FORCE SCALING   "
				<< consForceScaling << "\n" << endi;
     }
   }


   // external command forces

   if (extForcesOn) {
     iout << iINFO << "EXTERNAL COMMAND FORCES ACTIVE\n";
     iout << iINFO << "EXT FORCES COMMAND: " << extForcesCommand << "\n";
     iout << iINFO << "EXT COORD FILENAME: " << extCoordFilename << "\n";
     iout << iINFO << "EXT FORCE FILENAME: " << extForceFilename << "\n";
     iout << endi;
   }

    // QM command forces

    if (qmForcesOn) {
        iout << iINFO << "QM FORCES ACTIVE\n";
        if (qmParamPDBDefined){
            iout << iINFO << "QM PDB PARAMETER FILE: " << qmParamPDB << "\n";
        }
        iout << iINFO << "QM SOFTWARE: " << qmSoftware << "\n";

        if ( qmChrgMode == QMCHRGNONE )
            iout << iINFO << "QM ATOM CHARGES FROM QM SOFTWARE: NONE\n";
        if ( qmChrgMode == QMCHRGMULLIKEN )
            iout << iINFO << "QM ATOM CHARGES FROM QM SOFTWARE: MULLIKEN\n";
        if ( qmChrgMode == QMCHRGCHELPG )
            iout << iINFO << "QM ATOM CHARGES FROM QM SOFTWARE: CHELPG\n";

        iout << iINFO << "QM EXECUTABLE PATH: " << qmExecPath << "\n";
        iout << iINFO << "QM COLUMN: " << qmColumn << "\n";
        if (qmBondOn) {
            iout << iINFO << "QM BOND COLUMN: " << qmBondColumn << "\n";
            iout << iINFO << "QM WILL DETECT BONDS BETWEEN QM AND MM ATOMS.\n";
            if (qmBondDist) {
                iout << iINFO << "QM BOND COLUMN WILL DEFINE LINK AOTM DISTANCE.\n";
                if (qmBondValType == 1)
                    iout << iINFO << "QM BOND COLUMN HAS LENGTH INFORMATION.\n";
                else if (qmBondValType == 2)
                    iout << iINFO << "QM BOND COLUMN HAS RATIO INFORMATION.\n";
            }
            if (qmNoPC) {
                iout << iINFO << "MECHANICHAL EMBEDDING SELECTED."
                " BOND SCHEME WILL BE IGNORED!\n" << endi;
                qmBondScheme = QMSCHEMEZ1;
            }
            else {
                if (qmBondScheme == QMSCHEMECS)
                    iout << iINFO << "QM-MM BOND SCHEME: Charge Shift.\n";
                else if (qmBondScheme == QMSCHEMERCD)
                    iout << iINFO << "QM-MM BOND SCHEME: Redistributed Charge and Dipole.\n";
                else if (qmBondScheme == QMSCHEMEZ1)
                    iout << iINFO << "QM-MM BOND SCHEME: Z1.\n";
                else if (qmBondScheme == QMSCHEMEZ2)
                    iout << iINFO << "QM-MM BOND SCHEME: Z2.\n";
                else if (qmBondScheme == QMSCHEMEZ3)
                    iout << iINFO << "QM-MM BOND SCHEME: Z3.\n";
            }

        }

        if (qmChrgFromPSF) {
            iout << iINFO << "QM Will use PSF charges.\n";
        }

        iout << iINFO << "QM BASE DIRECTORY: " << qmBaseDir << "\n";

        if (qmPrepProcOn) {
            iout << iINFO << "QM PREPARATION PROCESS: " << qmPrepProc << "\n";
        }
        if (qmSecProcOn) {
            iout << iINFO << "QM SECONDARY PROCESS: " << qmSecProc << "\n";
        }

        current = config->find("QMConfigLine");
        for ( ; current; current = current->next ) {

            if ( strstr(current->data,"\n") ) {
                iout << iINFO << "QM configuration lines from NADM config file\n";
                continue;
            }

            iout << iINFO << "QM CONFIG LINE: " << current->data << "\n";

        }

        if (qmReplaceAll) {
            iout << iINFO << "QM FORCES WILL REPLACE ALL NAMD FORCES!\n";
        }

        if (qmNoPC)
            iout << iINFO << "QM NO POINT CHARGE: ON.\n";

        if (qmCustomPCSel)
            iout << iINFO << "QM CUSTOM POINT CHARGE SELECTION IS ACTIVATED\n";

        if (! qmNoPC && ! qmCustomPCSel)
            iout << iINFO << "QM POINT CHARGES WILL BE SELECTED EVERY "
            << qmPCSelFreq << " STEPS.\n";

        if (qmPCSwitchOn) {
            iout << iINFO << "QM Point Charge Switching: ON.\n";

            if (qmPCScheme == 1)
                iout << iINFO << "QM Point Charge SCHEME: none.\n";
            else if (qmPCScheme == 2)
                iout << iINFO << "QM Point Charge SCHEME: round.\n";
            else if (qmPCScheme == 3)
                iout << iINFO << "QM Point Charge SCHEME: zero.\n";
        }

        if (qmLSSOn) {
            iout << iINFO << "QM LIVE SOLVENT SELECTION IS ACTIVE.\n" ;
            iout << iINFO << "QM LIVE SOLVENT SELECTION FREQUENCY: "
            << qmLSSFreq << "\n" << endi;

            current = config->find("QMLSSSize");
            for ( ; current; current = current->next ) {
                iout << iINFO << "QM LIVE SOLVENT SELECTION SIZE (\"qmGrpID numMolecules\"): " << current->data << "\n";
            }

            if (! opts.defined("QMLWSResname"))
                strcpy(qmLSSResname,"TIP3");
            iout << iINFO << "QM LIVE SOLVENT SELECTION WILL USE RESIDUE TYPE: " << qmLSSResname << "\n" << endi;
        }

        iout << iINFO << "QM executions per node: " << qmSimsPerNode << "\n";

        iout << endi;
    }


   // gbis gbobc implicit solvent parameters

  if (GBISserOn) {
      GBISOn = 0;//turning gbis-ser on turns gbis-parallel off
     iout << iINFO<< "GBIS GENERALIZED BORN IMPLICIT SOLVENT ACTIVE (SERIAL)\n";
  }
  if (GBISOn) {
     iout << iINFO << "GBIS GENERALIZED BORN IMPLICIT SOLVENT ACTIVE\n";
  }
  if (GBISOn || GBISserOn) {
     iout << iINFO << "GBIS SOLVENT DIELECTRIC: " << solvent_dielectric<< "\n";
     iout << iINFO << "GBIS PROTEIN DIELECTRIC: " << dielectric<< "\n";
     iout <<iINFO<<"GBIS COULOMB RADIUS OFFSET: "<< coulomb_radius_offset<<" Ang\n";
     iout << iINFO << "GBIS ION CONCENTRATION: " << ion_concentration << " M\n";
     iout << iINFO << "GBIS DEBYE SCREENING LENGTH: " << 1.0/kappa << " Ang\n";
     iout << iINFO << "GBIS DELTA: " << gbis_delta << "\n";
     iout << iINFO << "GBIS BETA: " << gbis_beta << "\n";
     iout << iINFO << "GBIS GAMMA: " << gbis_gamma << "\n";
     iout << iINFO << "GBIS BORN RADIUS CUTOFF: " << alpha_cutoff << " Ang\n";
     iout << iINFO << "GBIS MAX BORN RADIUS: " << alpha_max << " Ang\n";
     iout << endi;
  }

  if (LCPOOn) {
    iout << iINFO << "SASA SURFACE TENSION: " << surface_tension<< " kcal/mol/Ang^2\n";
  }

   tclBCScript = 0;
   if (tclBCOn) {
     iout << iINFO << "TCL BOUNDARY FORCES ACTIVE\n";
     current = config->find("tclBCScript");
     if ( current ) {
       tclBCScript = current->data;
       iout << iINFO << "TCL BOUNDARY FORCES SCRIPT   " << current->data << "\n";
     }
       iout << iINFO << "TCL BOUNDARY FORCES ARGS     " << tclBCArgs << "\n";
     iout << endi;
   }

   // Global forces configuration

   globalForcesOn = ( tclForcesOn || freeEnergyOn || miscForcesOn ||
                      (IMDon && ! (IMDignore || IMDignoreForces)) || SMDOn || TMDOn ||
                      colvarsOn || symmetryOn || qmForcesOn );


   if (tclForcesOn)
   {
     iout << iINFO << "TCL GLOBAL FORCES ACTIVE\n";

     current = config->find("tclForcesScript");

     for ( ; current; current = current->next ) {

     if ( strstr(current->data,"\n") ) {
       iout << iINFO << "TCL GLOBAL FORCES SCRIPT INLINED IN CONFIG FILE\n";
       continue;
     }

     iout << iINFO << "TCL GLOBAL FORCES SCRIPT   " << current->data << "\n";

     }
     iout << endi;
   }

   if (miscForcesOn)
   {
     iout << iINFO << "MISC FORCES ACTIVE\n";

     current = config->find("miscForcesScript");

     for ( ; current; current = current->next ) {

     if ( strstr(current->data,"\n") ) {
       iout << iINFO << "MISC FORCES SCRIPT INLINED IN CONFIG FILE\n";
       continue;
     }

     iout << iINFO << "MISC FORCES SCRIPT   " << current->data << "\n";

     }
     iout << endi;
   }

   if (freeEnergyOn)
   {
     iout << iINFO << "FREE ENERGY PERTURBATION ACTIVE\n";

     current = config->find("freeEnergyConfig");

     for ( ; current; current = current->next ) {

     if ( strstr(current->data,"\n") ) {
       iout << iINFO << "FREE ENERGY PERTURBATION SCRIPT INLINED IN CONFIG FILE\n";
       continue;
     }

     iout << iINFO << "FREE ENERGY PERTURBATION SCRIPT   " << current->data << "\n";

     }
     iout << endi;
   }

   if (colvarsOn)
   {
     iout << iINFO << "COLLECTIVE VARIABLES CALCULATION REQUESTED\n";

     current = config->find ("colvarsConfig");
     for ( ; current; current = current->next ) {
       if ( strstr(current->data,"\n") ) {
         iout << iINFO << "COLLECTIVE VARIABLES CONFIGURATION INLINED IN CONFIG FILE\n";
         continue;
       }
       iout << iINFO << "COLLECTIVE VARIABLES CONFIGURATION   " << current->data << "\n";
     }

     current = config->find ("colvarsInput");
     for ( ; current; current = current->next ) {
       if ( strstr(current->data,"\n") ) {
         iout << iINFO << "COLLECTIVE VARIABLES RESTART INFORMATION INLINED IN CONFIG FILE\n";
         continue;
       }
       iout << iINFO << "COLLECTIVE VARIABLES RESTART INFORMATION   " << current->data << "\n";
     }

     iout << endi;
   }

   if (IMDon)
   {
     iout << iINFO << "INTERACTIVE MD ACTIVE\n";
     iout << iINFO << "INTERACTIVE MD PORT    " << IMDport << "\n";
     iout << iINFO << "INTERACTIVE MD FREQ    " << IMDfreq << "\n";
     if (IMDignore) {
        iout << iINFO << "INTERACTIVE MD WILL NOT INFLUENCE SIMULATION\n";
     } else {
       if (IMDignoreForces)
         {
            iout << iINFO << "INTERACTIVE FORCES ARE DISABLED\n";
            iout << iINFO << "PAUSE, RESUME, DETACH AND FINISH INTERACTIVE MD ARE ENABLED\n";
         }
       if (IMDwait) iout << iINFO << "WILL AWAIT INTERACTIVE MD CONNECTION\n";
     }
     iout << endi;
   }

   if (globalOn && !dihedralOn)
   {
      iout << iINFO << "GLOBAL INTEGRATION TEST MODE ACTIVE\n";
   }


   if (dihedralOn)
   {
      iout << iINFO << "DIHEDRAL ANGLE DYNAMICS ACTIVE\n";
      if (!COLDOn)
      {
         iout << iINFO << "*** DIHEDRAL ANGLE DYNAMICS IS HIGHLY EXPERIMENTAL ***\n";
         iout << iINFO << "PLEASE CONSIDER USING THE COLD OPTION AS WELL\n";
      }
   }

   // This function is so long that it exceeds the emacs brace
   // matching default max length of 25600 (a bit before this comment). We
   // should take this is a not too subtle hint about scale of this
   // violation of good coding practices.  Especially considering the
   // fact that this function still has about a thousand lines to go
   // before its done, and is doomed to grow with new features.

   if (COLDOn)
   {
      iout << iINFO << "COLD (CONSTRAINED OVERDAMPED LANGEVIN DYNAMICS) ACTIVE\n";

      iout << iINFO << "COLD TARGET TEMP       "
         << COLDTemp << "\n";

      iout << iINFO << "COLD COLLISION RATE    "
         << COLDRate << "\n";
   }

   if (cylindricalBCOn)
   {
    iout << iINFO << "CYLINDRICAL BOUNDARY CONDITIONS ACTIVE\n";
    iout << iINFO << "AXIS                     " << cylindricalBCAxis << "\n";
    iout << iINFO << "RADIUS #1                " << cylindricalBCr1 << "\n";
    iout << iINFO << "FORCE CONSTANT #1        " << cylindricalBCk1 << "\n";
    iout << iINFO << "EXPONENT #1              " << cylindricalBCexp1 << "\n";
    iout << iINFO << "LENGTH #1                " << cylindricalBCl1 << "\n";
    if (cylindricalBCr2 > 0.0)
    {
     iout << iINFO << "RADIUS #2               " << cylindricalBCr2 << "\n";
     iout << iINFO << "FORCE CONSTANT #2       " << cylindricalBCk2 << "\n";
     iout << iINFO << "EXPONENT #2             " << cylindricalBCexp2 << "\n";
     iout << iINFO << "LENGTH #2               " << cylindricalBCl2 << "\n";
    }
    iout << iINFO << "CYLINDER BOUNDARY CENTER(" << cylindricalCenter.x << ", "
             << cylindricalCenter.y << ", " << cylindricalCenter.z << ")\n";
    iout << endi;
  }

   if (sphericalBCOn)
   {
      iout << iINFO << "SPHERICAL BOUNDARY CONDITIONS ACTIVE\n";

      iout << iINFO << "RADIUS #1              "
         << sphericalBCr1 << "\n";
      iout << iINFO << "FORCE CONSTANT #1      "
         << sphericalBCk1 << "\n";
      iout << iINFO << "EXPONENT #1            "
         << sphericalBCexp1 << "\n";

      if (sphericalBCr2 > 0)
      {
        iout << iINFO << "RADIUS #2              "
              << sphericalBCr2 << "\n";
        iout << iINFO << "FORCE CONSTANT #2      "
            << sphericalBCk2 << "\n";
        iout << iINFO << "EXPONENT #2            "
        << sphericalBCexp2 << "\n";
      }

      iout << iINFO << "SPHERE BOUNDARY CENTER(" << sphericalCenter.x << ", "
               << sphericalCenter.y << ", " << sphericalCenter.z << ")\n";
      iout << endi;
   }

   if (eFieldOn)
   {
      iout << iINFO << "ELECTRIC FIELD ACTIVE\n";

      iout << iINFO << "E-FIELD VECTOR         ("
         << eField.x << ", " << eField.y
         << ", " << eField.z << ")\n";
      if ( eFieldNormalized ) iout << iINFO << "E-FIELD VECTOR IS SCALED BY CELL BASIS VECTORS\n";
      iout << iINFO << "E-FIELD FREQUENCY IS (1/ps) " << eFieldFreq << "\n";
      iout << iINFO << "E-FIELD PHASE IS     (deg)  " << eFieldPhase << "\n";

      iout << endi;
   }

      if (stirOn)
   {
      iout << iINFO << "STIRRING TORQUES ACTIVE\n";

      iout << iINFO << "STIR STARTING THETA   (deg)  "<< stirStartingTheta << "\n";
      iout << iINFO << "STIR ANGULAR VELOCITY (deg/ts)   " << stirVel <<"\n";
      iout << iINFO << "STIR FORCE HARMONIC SPRING CONSTANT "<< stirK << "\n";
      iout << iINFO << "STIR AXIS OF ROTATION (DIRECTION)      ("
         << stirAxis.x << ", " << stirAxis.y
         << ", " << stirAxis.z << ")\n";
	          iout << iINFO << "STIR PIVOT POINT (COORDINATE)           ("
         << stirPivot.x << ", " << stirPivot.y
         << ", " << stirPivot.z << ")\n";
      current = config->find("stirFilename");

      iout << iINFO << "STIR ATOMS AND ORIGINAL POSITIONS FROM FILE    " <<current ->data << '\n';
      current = config->find("stirredAtomsCol");
      iout << iINFO <<"STIR FILE COLUMN " << current ->data << '\n';
      iout << endi;
   }

   if (drudeOn)
   {
      iout << iINFO << "DRUDE MODEL DUAL THERMOSTAT IS ACTIVE\n";
      iout << iINFO << "DRUDE BOND TEMPERATURE " << drudeTemp << "\n";
      if (drudeDamping > 0.0) {
        iout << iINFO << "DRUDE DAMPING COEFFICIENT IS "
             << drudeDamping << " INVERSE PS\n";
      }
      if (drudeHardWallOn) {
        iout << iINFO << "DRUDE HARD WALL RESTRAINT IS ACTIVE FOR DRUDE BONDS\n";
        iout << iINFO << "DRUDE MAXIMUM BOND LENGTH BEFORE RESTRAINT IS   "
             << drudeBondLen << "\n";
      } else if (drudeBondConst > 0.0) {
        iout << iINFO << "DRUDE QUARTIC RESTRAINT IS ACTIVE FOR DRUDE BONDS\n";
        iout << iINFO << "DRUDE MAXIMUM BOND LENGTH BEFORE RESTRAINT IS   "
             << drudeBondLen << "\n";
        iout << iINFO << "DRUDE BOND RESTRAINT CONSTANT IS                "
             << drudeBondConst << "\n";
      }
      if (drudeNbtholeCut > 0.0) {
        iout << iINFO << "DRUDE NBTHOLE IS ACTIVE\n";
        iout << iINFO << "DRUDE NBTHOLE RADIUS IS   "
             << drudeNbtholeCut << "\n";
      }
   }

   if (langevinOn)
   {
      iout << iINFO << "LANGEVIN DYNAMICS ACTIVE\n";
      iout << iINFO << "LANGEVIN TEMPERATURE   "
         << langevinTemp << "\n";
      if (! langevin_useBAOAB) iout << iINFO << "LANGEVIN USING BBK INTEGRATOR\n";
      else  iout << iINFO << "LANGEVIN USING BAOAB INTEGRATOR\n"; // [!!] Info file
      if (langevinDamping > 0.0) {
	iout << iINFO << "LANGEVIN DAMPING COEFFICIENT IS "
		<< langevinDamping << " INVERSE PS\n";
	if (langevinHydrogen)
		iout << iINFO << "LANGEVIN DYNAMICS APPLIED TO HYDROGENS\n";
	else
		iout << iINFO << "LANGEVIN DYNAMICS NOT APPLIED TO HYDROGENS\n";
      } else {
	iout << iINFO << "LANGEVIN DAMPING COEFFICIENTS DETERMINED FROM FILES\n";
        current = config->find("langevinFile");
	if ( current ) iout << iINFO << "LANGEVIN DAMPING FILE:  " <<
          current->data << "\n";
        else iout << iINFO << "LANGEVIN DAMPING FILE IS COORDINATE PDB\n";
        current = config->find("langevinCol");
	if ( current ) iout << iINFO << "LANGEVIN DAMPING COLUMN:  " <<
          current->data << "\n";
        else iout << iINFO << "LANGEVIN DAMPING COLUMN:  DEFAULT (4TH, O)\n";
      }
      iout << endi;
   }

   // BEGIN LA
   if (loweAndersenOn)
   {
      iout << iINFO << "LOWE-ANDERSEN DYNAMICS ACTIVE\n";
      iout << iINFO << "LOWE-ANDERSEN TEMPERATURE     "
         << loweAndersenTemp << " K\n";
      iout << iINFO << "LOWE-ANDERSEN RATE            "
         << loweAndersenRate << " INVERSE PS\n";
      iout << iINFO << "LOWE-ANDERSEN CUTOFF          "
         << loweAndersenCutoff << " ANGSTROMS\n";
      iout << endi;
   }
   // END LA

   if (tCoupleOn)
   {
      iout << iINFO << "TEMPERATURE COUPLING ACTIVE\n";
      iout << iINFO << "COUPLING TEMPERATURE   "
         << tCoupleTemp << "\n";
      iout << endi;
   }

   if (stochRescaleOn)
   {
      iout << iINFO << "STOCHASTIC RESCALING ACTIVE\n";
      iout << iINFO << "STOCHASTIC RESCALING TEMPERATURE "
           << stochRescaleTemp << " K\n";
      iout << iINFO << "STOCHASTIC RESCALING PERIOD "
           << stochRescalePeriod << " PS\n";
      iout << iINFO << "STOCHASTIC RESCALING WILL OCCUR EVERY "
           << stochRescaleFreq << " STEPS\n";
      iout << endi;
   }

   if (minimizeOn)
   {
      iout << iINFO << "OLD STYLE MINIMIZATION ACTIVE\n";
      iout << endi;
   }

   if (minimizeCGOn)
   {
      iout << iINFO << "CONJUGATE GRADIENT MINIMIZATION ACTIVE\n";
      iout << iINFO << "LINE MINIMIZATION GOAL = " << minLineGoal << "\n";
      iout << iINFO << "BABY STEP SIZE = " << minBabyStep << "\n";
      iout << iINFO << "TINY STEP SIZE = " << minTinyStep << "\n";
      iout << endi;
   }

   if (maximumMove)
   {
      iout << iINFO << "MAXIMUM MOVEMENT       "
         << maximumMove << "\n";
      iout << endi;
   }

   if (rescaleFreq > 0)
   {
     iout << iINFO << "VELOCITY RESCALE FREQ  "
        << rescaleFreq << "\n";
     iout << iINFO << "VELOCITY RESCALE TEMP  "
        << rescaleTemp << "\n";
     iout << endi;
   }

   if (reassignFreq > 0)
   {
     iout << iINFO << "VELOCITY REASSIGNMENT FREQ  "
        << reassignFreq << "\n";
     iout << iINFO << "VELOCITY REASSIGNMENT TEMP  "
        << reassignTemp << "\n";
     if ( reassignIncr != 0. )
       iout << iINFO << "VELOCITY REASSIGNMENT INCR  "
        << reassignIncr << "\n";
     if ( reassignHold != 0. )
       iout << iINFO << "VELOCITY REASSIGNMENT HOLD  "
        << reassignHold << "\n";
     iout << endi;
   }

   if ((int)berendsenPressureOn + (int)langevinPistonOn + (int)multigratorOn > 1)
   {
      NAMD_die("Multiple pressure control algorithms selected!\n");
   }

   if (excludeFromPressure) {
     iout << iINFO << "EXCLUDE FROM PRESSURE ACTIVE\n";
   }
   if (useConstantArea && useConstantRatio) {
     NAMD_die("useConstantArea and useConstantRatio are mutually exclusive.\n");
   }
   if (useConstantRatio && !useFlexibleCell) {
     NAMD_die("useConstantRatio requires useFlexibleCell.\n");
   }
   if (useConstantArea && surfaceTensionTarget) {
     NAMD_die("surfaceTensionTarget and useConstantArea are mutually exclusive.\n");
   }
   if (useConstantArea && !useFlexibleCell) {
     NAMD_die("useConstantArea requires useFlexibleCell.\n");
   }

   if (berendsenPressureOn || langevinPistonOn) {
     if (rigidBonds != RIGID_NONE && useGroupPressure == FALSE) {
       useGroupPressure = TRUE;
       iout << iWARN << "Option useGroupPressure is being enabled "
            << "due to pressure control with rigidBonds.\n" << endi;
     }
   }

   if (berendsenPressureOn)
   {
     if ( ! opts.defined("BerendsenPressureFreq") ) {
	berendsenPressureFreq = nonbondedFrequency;
	if ( fullElectFrequency )
		berendsenPressureFreq = fullElectFrequency;
     }
     if ( (berendsenPressureFreq % nonbondedFrequency) || ( fullElectFrequency
		&& (berendsenPressureFreq % fullElectFrequency) ) )
	NAMD_die("berendsenPressureFreq must be a multiple of both fullElectFrequency and nonbondedFrequency\n");
     iout << iINFO << "BERENDSEN PRESSURE COUPLING ACTIVE\n";
     iout << iINFO << "    TARGET PRESSURE IS "
        << berendsenPressureTarget << " BAR\n";
     iout << iINFO << "    COMPRESSIBILITY ESTIMATE IS "
        << berendsenPressureCompressibility << " BAR^(-1)\n";
     iout << iINFO << "    RELAXATION TIME IS "
        << berendsenPressureRelaxationTime << " FS\n";
     iout << iINFO << "    APPLIED EVERY "
        << berendsenPressureFreq << " STEPS\n";
     iout << iINFO << "    PRESSURE CONTROL IS "
	<< (useGroupPressure?"GROUP":"ATOM") << "-BASED\n";
     iout << endi;
     berendsenPressureTarget /= PRESSUREFACTOR;
     berendsenPressureCompressibility *= PRESSUREFACTOR;
   }

   if (langevinPistonOn)
   {
     iout << iINFO << "LANGEVIN PISTON PRESSURE CONTROL ACTIVE\n";
     iout << iINFO << "       TARGET PRESSURE IS "
        << langevinPistonTarget << " BAR\n";
     iout << iINFO << "    OSCILLATION PERIOD IS "
        << langevinPistonPeriod << " FS\n";
     iout << iINFO << "            DECAY TIME IS "
        << langevinPistonDecay << " FS\n";
     iout << iINFO << "    PISTON TEMPERATURE IS "
        << langevinPistonTemp << " K\n";
     iout << iINFO << "      PRESSURE CONTROL IS "
	<< (useGroupPressure?"GROUP":"ATOM") << "-BASED\n";
     iout << iINFO << "   INITIAL STRAIN RATE IS "
        << strainRate << "\n";
     iout << endi;
     langevinPistonTarget /= PRESSUREFACTOR;
   }

    if (multigratorOn) {
      multigratorPressureTarget /= PRESSUREFACTOR;
    }

   if (berendsenPressureOn || langevinPistonOn) {
     iout << iINFO << "      CELL FLUCTUATION IS "
	    << (useFlexibleCell?"AN":"") << "ISOTROPIC\n";
     if (useConstantRatio)
       iout << iINFO << "    SHAPE OF CELL IS CONSTRAINED IN X-Y PLANE\n";
     if (useConstantArea)
       iout << iINFO << "    CONSTANT AREA PRESSURE CONTROL ACTIVE\n";
   }

   if (surfaceTensionTarget != 0)
   {
     iout << iINFO << "SURFACE TENSION CONTROL ACTIVE\n";
     iout << iINFO << "      TARGET SURFACE TENSION IS "
          << surfaceTensionTarget << " DYN/CM\n";
     iout << endi;
     // multiply by 100 to convert from dyn/cm to bar-Angstroms, then divide
     // by PRESSURE factor to convert bar to NAMD internal pressure units.
     surfaceTensionTarget *= 100.0 / PRESSUREFACTOR;
   }

   if (pressureProfileOn) {
     if ((berendsenPressureOn || langevinPistonOn) && !dcdUnitCell) {
#if 1
       iout << iWARN << "Turning on dcdUnitCell so that trajectory files contain unit cell data.\n" << endi;
       dcdUnitCell = 1;
#else
       NAMD_die("Sorry, pressure profile not implemented for constant pressure.");
#endif
     }
     // if Ewald is on, only calculate Ewald
     if (pressureProfileEwaldOn)
       pressureProfileOn = 0;

     if (pressureProfileSlabs < 1)
       NAMD_die("pressureProfileSlabs must be positive.");
     iout << iINFO << "PRESSURE PROFILE CALCULATIONS ACTIVE\n";
     iout << iINFO << "      NUMBER OF SLABS: " << pressureProfileSlabs << "\n";
     iout << iINFO << "      SLAB THICKNESS: " << cellBasisVector3.z / pressureProfileSlabs
                   << "\n";
     iout << iINFO << "      TIMESTEPS BETWEEN DATA OUTPUT: "
                   << pressureProfileFreq << "\n";
     iout << iINFO << "      NUMBER OF ATOM TYPES: " << pressureProfileAtomTypes << "\n";
     iout << endi;
   } else {
     pressureProfileEwaldOn = 0;
     pressureProfileAtomTypes = 1;
   }

   if (accelMDOn) {
     iout << iINFO << "ACCELERATED MD ACTIVE\n";

     if ( accelMDdual) {
        accelMDdihe = FALSE;
        iout << iINFO << "APPLYING DUAL BOOST\n";
     }
     else if ( accelMDdihe ) {
        iout << iINFO << "BOOSTING DIHEDRAL POTENTIAL\n";
     } else {
        iout << iINFO << "BOOSTING TOTAL POTENTIAL\n";
     }

     if(accelMDG){
	 switch(accelMDGiE) {
	     case 1:
		 iout << iINFO << "accelMDG THRESHOLD ENERGY SET TO LOWER BOUND Vmax\n";
		 break;
	     case 2:
		 iout << iINFO << "accelMDG THRESHOLD ENERGY SET TO UPPER BOUND Vmin+(Vmax-Vmin)/k0\n";
		 break;
	 }
	 if(accelMDGRestart)
	     iout << iINFO << "accelMDG USING RESTART FILE " << accelMDGRestartFile << "\n";
	 if(accelMDGresetVaftercmd)
	     iout << iINFO << "accelMDG WILL RESET STATISTICS AFTER FIRST CMD STEPS\n";

	 iout << iINFO << "accelMDG " << accelMDGcMDSteps << " CONVENTIONAL MD STEPS "
	     << "(WITH " << accelMDGcMDPrepSteps << " PREPARATION STEPS)\n";
	 if(accelMDGcMDSteps == 0)
		 iout << iINFO << "(accelMDGcMDPrepSteps is set to zero automatically)\n";

	 iout << iINFO << "accelMDG " << accelMDGEquiSteps << " EQUILIBRATION STEPS "
	     << "(WITH " << accelMDGEquiPrepSteps << " PREPARATION STEPS)\n";
	 if(accelMDGEquiSteps == 0)
		 iout << iINFO << "(accelMDGEquiPrepSteps is set to zero automatically)\n";

         if(accelMDGStatWindow > 0)
             iout << iINFO << "accelMDG WILL RESET AVERAGE AND STANDARD DEVIATION EVERY " << accelMDGEquiSteps << " STEPS\n";
         else
             iout << iINFO << "accelMDG WILL NOT RESET AVERAGE AND STANDARD DEVIATION\n";

	 if(accelMDdihe)
	     iout << iINFO << "accelMDGSigma0D: " << accelMDGSigma0D << " KCAL/MOL\n";
	 else if(accelMDdual)
	     iout << iINFO << "accelMDGSigma0P: " << accelMDGSigma0P << " KCAL/MOL, "
		 << "accelMDGSigma0D: " << accelMDGSigma0D << " KCAL/MOL\n";
	 else
	     iout << iINFO << "accelMDGSigma0P: " << accelMDGSigma0P << " KCAL/MOL\n";
     }
     else{
	 iout << iINFO << "accelMDE: " << accelMDE << " KCAL/MOL, accelMDalpha: " << accelMDalpha << " KCAL/MOL\n";
	 if (accelMDdual) {
	     iout << iINFO << "accelMDTE: " << accelMDTE << " KCAL/MOL, "
		 << "accelMDTalpha: " << accelMDTalpha << " KCAL/MOL\n";
	 }
     }
     if ( accelMDLastStep > 0) {
        iout << iINFO << "accelMD WILL BE DONE FROM STEP " << accelMDFirstStep << " TO STEP " << accelMDLastStep << "\n";
     } else {
        iout << iINFO << "accelMD WILL BE DONE FROM STEP " << accelMDFirstStep << " TO THE END OF THE SIMULATION \n";
     }
     iout << iINFO << "accelMD OUTPUT FREQUENCY " << accelMDOutFreq << "\n";
     iout << endi;
   }

   if (adaptTempOn) {
     iout << iINFO << "ADAPTIVE TEMPERING ACTIVE:\n";
     iout << iINFO << "      OUTPUT FREQUENCY: " << adaptTempOutFreq << "\n";
     iout << iINFO << "      TEMPERATURE UPDATE FREQUENCY: " << adaptTempFreq << "\n";
     if ( adaptTempLastStep > 0 )
        iout << iINFO << "      ADAPTIVE TEMPERING WILL BE DONE FROM STEP " << adaptTempFirstStep  << " TO " << adaptTempLastStep << "\n";
     else
        iout << iINFO << "      ADAPTIVE TEMPERING WILL BE DONE FROM STEP " << adaptTempFirstStep << "\n";
     if ( adaptTempLangevin )
        iout << iINFO << "      ADAPTIVE TEMPERING COUPLED TO LANGEVIN THERMOSTAT\n";
     if ( adaptTempRescale )
        iout << iINFO << "      ADAPTIVE TEMPERING COUPLED TO VELOCITY RESCALING\n";
     if (adaptTempRestartFreq > 0) {
        iout << iINFO << "      WRITING RESTART INFORMATION TO " << adaptTempRestartFile << " EVERY " << adaptTempRestartFreq << " STEPS\n";
     }

   }

   if (FMAOn)
   {
     iout << iINFO << "FMA ACTIVE\n";
     iout << iINFO << "FMA THETA              "
        << fmaTheta << "\n";
     iout << endi;
   }

   FFTWWisdomString = 0;
   if (PMEOn)
   {
     iout << iINFO << "PARTICLE MESH EWALD (PME) ACTIVE\n";
     iout << iINFO << "PME TOLERANCE               "
	<< PMETolerance << "\n";
     iout << iINFO << "PME EWALD COEFFICIENT       "
	<< PMEEwaldCoefficient << "\n";
     iout << iINFO << "PME INTERPOLATION ORDER     "
	<< PMEInterpOrder << "\n";
     iout << iINFO << "PME GRID DIMENSIONS         "
	<< PMEGridSizeX << " "
	<< PMEGridSizeY << " "
	<< PMEGridSizeZ << "\n";
     iout << iINFO << "PME MAXIMUM GRID SPACING    "
	<< PMEGridSpacing << "\n";
     if ( PMEBarrier ) {
       iout << iINFO << "PME BARRIER ENABLED\n";
     }
     if ( PMEOffload ) {
       iout << iINFO << "PME RECIPROCAL SUM OFFLOADED TO GPU\n";
     }
     iout << endi;
     if ( useDPME ) iout << iINFO << "USING OLD DPME CODE\n";
#ifdef NAMD_FFTW
     else if ( FFTWUseWisdom ) {  // handle FFTW wisdom
#ifdef NAMD_FFTW_3
       iout << iINFO << "Attempting to read FFTW data from system" <<"\n" <<endi;
       fftwf_import_system_wisdom();
#endif
       if (! opts.defined("FFTWWisdomFile")) {
         strcpy(FFTWWisdomFile,"FFTW_NAMD_");
         strcat(FFTWWisdomFile,NAMD_VERSION);
	 strcat(FFTWWisdomFile,"_");
	 strcat(FFTWWisdomFile,NAMD_PLATFORM);
#ifdef NAMD_FFTW_3
	 strcat(FFTWWisdomFile,"_FFTW3");
#endif
	 strcat(FFTWWisdomFile,".txt");
       }

       iout << iINFO << "Attempting to read FFTW data from "
		<< FFTWWisdomFile << "\n" << endi;
       FILE *wisdom_file = fopen(FFTWWisdomFile,"r");
       if ( wisdom_file ) {
#ifdef NAMD_FFTW_3
	 fftwf_import_wisdom_from_file(wisdom_file);
#else
	 fftw_import_wisdom_from_file(wisdom_file);
#endif
	 fclose(wisdom_file);
       }
       int nrp = 1;

       // rules based on work available
       int minslices = PMEMinSlices;
       int dimx = PMEGridSizeX;
       int nrpx = ( dimx + minslices - 1 ) / minslices;
       if ( nrpx > nrp ) nrp = nrpx;
       int dimy = PMEGridSizeY;
       int nrpy = ( dimy + minslices - 1 ) / minslices;
       if ( nrpy > nrp ) nrp = nrpy;

       // rules based on processors available
       int nrpp = CkNumPes();
       // if ( nrpp > 32 ) nrpp = 32;  // cap to limit messages
       if ( nrpp < nrp ) nrp = nrpp;

       // user override
       int nrps = PMEProcessors;
       if ( nrps > CkNumPes() ) nrps = CkNumPes();
       if ( nrps > 0 ) nrp = nrps;

       // make sure there aren't any totally empty processors
       int bx = ( dimx + nrp - 1 ) / nrp;
       int nrpbx = ( dimx + bx - 1 ) / bx;
       int by = ( dimy + nrp - 1 ) / nrp;
       int nrpby = ( dimy + by - 1 ) / by;
       nrp = ( nrpby > nrpbx ? nrpby : nrpbx );
       if ( bx != ( dimx + nrp - 1 ) / nrp )
         NAMD_bug("Error in selecting number of PME processors.");
       if ( by != ( dimy + nrp - 1 ) / nrp )
         NAMD_bug("Error in selecting number of PME processors.");

       // numGridPes = nrpbx;
       // numTransPes = nrpby;
       // numRecipPes = nrp;
       int block2 = (PMEGridSizeY + nrp - 1) / nrp;
       int block2_min = PMEGridSizeY % block2;
       if ( ! block2_min ) block2_min = block2;
       int dim3 = 2 * (PMEGridSizeZ/2 + 1);

       int n[3]; n[0] = PMEGridSizeX; n[1] = PMEGridSizeY; n[2] = PMEGridSizeZ;
       fftw_complex *work = new fftw_complex[n[0]];
       float *grid1 = (float *) fftwf_malloc(sizeof(float) *n[1]*dim3);
       float *grid2 = (float *) fftwf_malloc(sizeof(float) *n[0]*block2*dim3*2);
       iout << iINFO << "Optimizing 6 FFT steps.  1..." << endi;
#ifdef NAMD_FFTW_3
       int fftwFlags = FFTWPatient ? FFTW_PATIENT  : FFTWEstimate ? FFTW_ESTIMATE  : FFTW_MEASURE ;
       int planLineSizes[1];
       planLineSizes[0]=n[2];
       int nx= n[0];
       int ny=block2;
       int sizeLines=nx*ny;
       int zdim = dim3;
       int nz=zdim;
       int xStride=block2 * dim3 / 2;
       fftwf_destroy_plan(
			 fftwf_plan_many_dft_r2c(1, planLineSizes, sizeLines,
						 (float *) grid2, NULL, 1,
						 dim3,
						 (fftwf_complex *) grid2,
						 NULL, 1,
						 dim3/2,
						 fftwFlags));

       iout << " 2..." << endi;
       fftwf_destroy_plan(
			 fftwf_plan_many_dft_c2r(1, planLineSizes, sizeLines,
						 (fftwf_complex *) grid2,
						 NULL, 1,
						 dim3/2,
						 (float *) grid2, NULL, 1,
						 dim3,
						 fftwFlags));
       iout << " 3..." << endi;
       sizeLines=nz;
       planLineSizes[0]=block2;
       fftwf_destroy_plan(fftwf_plan_many_dft(1, planLineSizes, sizeLines,
					      (fftwf_complex *) grid2, NULL,
					      sizeLines, 1,
					      (fftwf_complex *) grid2, NULL,
					      sizeLines, 1,
					      FFTW_FORWARD,
					      fftwFlags));
       iout << " 4..." << endi;
       fftwf_destroy_plan(fftwf_plan_many_dft(1, planLineSizes, sizeLines,
					      (fftwf_complex *) grid2, NULL,
					      sizeLines, 1,
					      (fftwf_complex *) grid2, NULL,
					      sizeLines, 1,
					      FFTW_FORWARD,
					      fftwFlags));
       iout << " 5..." << endi;
       sizeLines=ny*nz;
       planLineSizes[0]=n[0];
       fftwf_destroy_plan(fftwf_plan_many_dft(1, planLineSizes, sizeLines,
					      (fftwf_complex *) grid2, NULL,
					      sizeLines, 1,
					      (fftwf_complex *) grid2, NULL,
					      sizeLines, 1,
					      FFTW_FORWARD,
					      fftwFlags));
       iout << " 6..." << endi;
       fftwf_destroy_plan(fftwf_plan_many_dft(1, planLineSizes, sizeLines,
					      (fftwf_complex *) grid2, NULL,
					      sizeLines, 1,
					      (fftwf_complex *) grid2, NULL,
					      sizeLines, 1,
					      FFTW_BACKWARD,
					      fftwFlags));

#else
       rfftwnd_destroy_plan( rfftwnd_create_plan_specific(
	 2, n+1, FFTW_REAL_TO_COMPLEX,
	 ( FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	 | FFTW_IN_PLACE | FFTW_USE_WISDOM, grid1, 1, 0, 0) );
       iout << " 2..." << endi;
       fftw_destroy_plan( fftw_create_plan_specific(n[0], FFTW_REAL_TO_COMPLEX,
	 ( FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	 | FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) grid2,
	 block2*dim3/2, work, 1) );
       iout << " 3..." << endi;
       fftw_destroy_plan( fftw_create_plan_specific(n[0], FFTW_REAL_TO_COMPLEX,
	 ( FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	 | FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) grid2,
	 block2_min*dim3/2, work, 1) );
       iout << " 4..." << endi;
       fftw_destroy_plan( fftw_create_plan_specific(n[0], FFTW_COMPLEX_TO_REAL,
	 ( FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	 | FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) grid2,
	 block2*dim3/2, work, 1) );
       iout << " 5..." << endi;
       fftw_destroy_plan( fftw_create_plan_specific(n[0], FFTW_COMPLEX_TO_REAL,
	 ( FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	 | FFTW_IN_PLACE | FFTW_USE_WISDOM, (fftw_complex *) grid2,
	 block2_min*dim3/2, work, 1) );
       iout << " 6..." << endi;
       rfftwnd_destroy_plan( rfftwnd_create_plan_specific(
	 2, n+1, FFTW_COMPLEX_TO_REAL,
	 ( FFTWEstimate ? FFTW_ESTIMATE : FFTW_MEASURE )
	 | FFTW_IN_PLACE | FFTW_USE_WISDOM, grid1, 1, 0, 0) );
#endif
       iout << "   Done.\n" << endi;
       delete [] work;
       fftwf_free(grid1);
       fftwf_free(grid2);

#ifdef NAMD_FFTW_3
       FFTWWisdomString = fftwf_export_wisdom_to_string();
#else
       FFTWWisdomString = fftw_export_wisdom_to_string();
#endif

      if ( FFTWWisdomString && (CmiNumPartitions() == 1) ) {
       iout << iINFO << "Writing FFTW data to "
		<< FFTWWisdomFile << "\n" << endi;
       wisdom_file = fopen(FFTWWisdomFile,"w");
       if ( wisdom_file ) {
#ifdef NAMD_FFTW_3
	 fftwf_export_wisdom_to_file(wisdom_file);
#else
	 fftw_export_wisdom_to_file(wisdom_file);
#endif
	 fclose(wisdom_file);
       }
      }
     }
#endif
     iout << endi;
   }
   if (fullDirectOn)
   {
     iout << iINFO << "DIRECT FULL ELECTROSTATIC CALCULATIONS ACTIVE\n";
     iout << endi;
   }

   // MSM configure
   if (MSMOn)
   {
     // check MSMQuality
     enum { LO=0, MEDLO, MED, MEDHI, HI };

     // MSMApprox
     enum { CUBIC=0, QUINTIC, QUINTIC2,
       SEPTIC, SEPTIC3, NONIC, NONIC4, C1HERMITE, NUM_APPROX };

     // MSMSplit
     enum { TAYLOR2=0, TAYLOR3, TAYLOR4,
       TAYLOR5, TAYLOR6, TAYLOR7, TAYLOR8, NUM_SPLIT };

     if (MSMApprox || MSMSplit) {  // take these definitions
       if (MSMApprox < 0 || MSMApprox >= NUM_APPROX) {
         NAMD_die("MSM: unknown approximation requested (MSMApprox)");
       }
       if (MSMSplit < 0 || MSMSplit >= NUM_SPLIT) {
         NAMD_die("MSM: unknown splitting requested (MSMSplit)");
       }
     }
     else {  // otherwise use MSMQuality to set MSMApprox and MSMSplit
       switch (MSMQuality) {
         case LO:
           MSMApprox = CUBIC;
           MSMSplit = TAYLOR2;
           break;
         case MEDLO:
           MSMApprox = C1HERMITE;
           MSMSplit = TAYLOR3;
           break;
         case MED:
           MSMApprox = QUINTIC;
           MSMSplit = TAYLOR3;
           break;
         case MEDHI:
           MSMApprox = SEPTIC;
           MSMSplit = TAYLOR4;
           break;
         case HI:
           MSMApprox = NONIC;
           MSMSplit = TAYLOR5;
           break;
         default:
           NAMD_die("MSM: unknown quality requested (MSMQuality)");
       }
     }

     iout << iINFO
       << "MULTILEVEL SUMMATION METHOD (MSM) FOR ELECTROSTATICS ACTIVE\n";
     if (MsmSerialOn) {
       iout << iINFO
         << "PERFORMING SERIAL MSM CALCULATION FOR LONG-RANGE PART\n";
     }
     const char *approx_str, *split_str;
     switch (MSMApprox) {
       case CUBIC:    approx_str = "C1 CUBIC";   break;
       case QUINTIC:  approx_str = "C1 QUINTIC"; break;
       case QUINTIC2: approx_str = "C2 QUINTIC"; break;
       case SEPTIC:   approx_str = "C1 SEPTIC";  break;
       case SEPTIC3:  approx_str = "C3 SEPTIC";  break;
       case NONIC:    approx_str = "C1 NONIC";   break;
       case NONIC4:   approx_str = "C4 NONIC";   break;
       case C1HERMITE:approx_str = "C1 HERMITE"; break;
       default:       approx_str = "UNKNOWN";    break;
     }
     switch (MSMSplit) {
       case TAYLOR2:  split_str = "C2 TAYLOR";   break;
       case TAYLOR3:  split_str = "C3 TAYLOR";   break;
       case TAYLOR4:  split_str = "C4 TAYLOR";   break;
       case TAYLOR5:  split_str = "C5 TAYLOR";   break;
       case TAYLOR6:  split_str = "C6 TAYLOR";   break;
       case TAYLOR7:  split_str = "C7 TAYLOR";   break;
       case TAYLOR8:  split_str = "C8 TAYLOR";   break;
       default:       split_str = "UNKNOWN";     break;
     }
     iout << iINFO
       << "MSM WITH " << approx_str << " INTERPOLATION "
       << "AND " << split_str << " SPLITTING\n"
       << endi;

   } // end MSM configure
   if (FMMOn)
   {
#ifdef FMM_SOLVER
     iout << iINFO << "FAST MULTIPOLE METHOD (FMM) FOR ELECTROSTATICS ACTIVE\n";
     iout << iINFO << "PERFORMING SERIAL FMM CALCULATION\n";
     iout << iINFO << "FMM LEVELS = " << FMMLevels << "\n";
     iout << iINFO << "FMM PADDING = " << FMMPadding << " ANGSTROMS\n";
     iout << endi;
#else
     NAMD_die("Must link to FMM library to use FMM\n");
#endif
   }

   if ( FMAOn || PMEOn || MSMOn || fullDirectOn || GBISOn || FMMOn )
   {
     iout << iINFO << "FULL ELECTROSTATIC EVALUATION FREQUENCY      "
	<< fullElectFrequency << "\n";
     iout << endi;

     if ( ( outputEnergies % fullElectFrequency ) &&
          ( fullElectFrequency % outputEnergies ) )
	NAMD_die("Either outputEnergies must be a multiple of fullElectFrequency or vice versa.\n");
   }

  if (MTSAlgorithm == NAIVE)
  {
    iout << iINFO << "USING NAIVE (CONSTANT FORCE) MTS SCHEME.\n" << endi;
  }
  if (MTSAlgorithm == VERLETI )
  {
    iout << iINFO << "USING VERLET I (r-RESPA) MTS SCHEME.\n" << endi;
  }

   if (longSplitting == SHARP)
  iout << iINFO << "SHARP SPLITTING OF LONG RANGE ELECTROSTATICS\n";
   else if (longSplitting == XPLOR)
  iout << iINFO << "XPLOR SPLITTING OF LONG RANGE ELECTROSTATICS\n";
   else if (longSplitting == C1)
  iout << iINFO << "C1 SPLITTING OF LONG RANGE ELECTROSTATICS\n";
   else if (longSplitting == C2)
  iout << iINFO << "C2 SPLITTING OF LONG RANGE ELECTROSTATICS\n";

   if (splitPatch == SPLIT_PATCH_POSITION)
  iout << iINFO << "PLACING ATOMS IN PATCHES BY POSITION\n";
   else if (splitPatch == SPLIT_PATCH_HYDROGEN)
  iout << iINFO << "PLACING ATOMS IN PATCHES BY HYDROGEN GROUPS\n";

   iout << endi;

   if (mollyOn)
   {
     iout << iINFO << "SLOW FORCE MOLLIFICATION : \n";
     iout << iINFO << "         ERROR TOLERANCE : " << mollyTol << "\n";
     iout << iINFO << "          MAX ITERATIONS : " << mollyIter << "\n";
     iout << endi;
   }

   if (rigidBonds != RIGID_NONE)
   {
     iout << iINFO << "RIGID BONDS TO HYDROGEN : ";
     if (rigidBonds == RIGID_ALL)    iout << "ALL\n";
     if (rigidBonds == RIGID_WATER)  iout << "WATER\n";
     iout << iINFO << "        ERROR TOLERANCE : " << rigidTol << "\n";
     iout << iINFO << "         MAX ITERATIONS : " << rigidIter << "\n";
     if (useSettle) iout << iINFO << "RIGID WATER USING SETTLE ALGORITHM\n";
     iout << endi;
   }


   if (nonbondedFrequency != 1)
   {
     iout << iINFO << "NONBONDED FORCES EVALUATED EVERY " << nonbondedFrequency << " STEPS\n";
   }

   iout << iINFO << "RANDOM NUMBER SEED     "
      << randomSeed << "\n";

   iout << endi;

   iout << iINFO << "USE HYDROGEN BONDS?    ";
   if (HydrogenBonds)
   {
  iout << "YES\n" << endi;
  iout << iINFO << "USE ANTECEDENT ATOMS?  ";
  iout << (useAntecedent ? "YES" : "NO");
        iout << "\nHB DIST CUT, ON, OFF   ";
  iout << daCutoffDist << " , " << daOnDist << " , " << daOffDist;
        iout << "\nHB ANGLE CUT, ON, OFF  ";
  iout << dhaCutoffAngle << " , " << dhaOnAngle << " , ";
  iout << dhaOffAngle;
        iout << "\nHB ATT, REP exponents  ";
  iout << distAttExp << " , " << distRepExp;
        iout << "\nHB AA, HA exponents    ";
  iout << aaAngleExp << " , " << haAngleExp;
  iout << "\n" << endi;
   }
   else
   {
  iout << "NO\n" << endi;
   }

// If this is AMBER, then print AMBER options

   if (amberOn)
   { iout << iINFO << "Using AMBER format force field!\n";
     current = config->find("parmfile");
     iout << iINFO << "AMBER PARM FILE        " << current->data << '\n';
     if (opts.defined("coordinates"))
     { current = config->find("coordinates");
       iout << iINFO << "COORDINATE PDB         " << current->data << '\n';
     }
     else
     { current = config->find("ambercoor");
       iout << iINFO << "AMBER COORDINATE FILE  " << current->data << '\n';
     }
     if (readExclusions)
       iout << iINFO << "Exclusions will be read from PARM file!\n";
     else
       iout << iINFO << "Exclusions in PARM file will be ignored!\n";
     iout << iINFO << "SCNB (VDW SCALING)     " << vdwscale14 << "\n" << endi;
   }
   else if(gromacsOn)
   {
     iout << iINFO << "Using GROMACS format force field!\n";

     current = config->find("grotopfile");
     // it should be defined, but, just in case...
     if (current == NULL)
       NAMD_die("no GROMACS topology file defined!?");
     iout << iINFO << "GROMACS TOPO FILE        " << current->data << '\n';

     // XXX handle the two types of coordinates more gracefully
     current = config->find("grocoorfile");
     if (current == NULL) {
       current = config->find("coordinates");
       if (current == NULL) {
	 NAMD_die("no coordinate file defined!?");
       }
     }
     iout << iINFO << "GROMACS COOR FILE        " << current->data << '\n'
	  << endi;

   }
   else {
     if ( !usePluginIO ) {
       if ( current = config->find("coordinates") )
       iout << iINFO << "COORDINATE PDB         " << current->data << '\n' << endi;
     }

     current = config->find("structure");

     iout << iINFO << "STRUCTURE FILE         "
        << current->data << "\n" << endi;

     if (cosAngles)
     {
       iout << iINFO << "COSANGLES ON. SOME ANGLES WILL BE COSINE-BASED\n" << endi;
     }

     //****** BEGIN CHARMM/XPLOR type changes
     if (paraTypeXplorOn)
     {
       iout << iINFO << "PARAMETER file: XPLOR format! (default) \n" << endi;
     }
     else if (paraTypeCharmmOn)
     {
       iout << iINFO << "PARAMETER file: CHARMM format! \n" << endi;
     }
     //****** END CHARMM/XPLOR type changes

     current = config->find("parameters");

     while (current != NULL)
     {
       iout << iINFO << "PARAMETERS             "
          << current->data << "\n" << endi;
       current = current->next;
     }
   }

     iout << iINFO << "USING " <<
        ( vdwGeometricSigma ? "GEOMETRIC" : "ARITHMETIC" ) <<
        " MEAN TO COMBINE L-J SIGMA PARAMETERS\n" << endi;

   if (opts.defined("bincoordinates"))
   {
     current = config->find("bincoordinates");

     iout << iINFO << "BINARY COORDINATES     "
              << current->data << "\n";
   }

#ifdef MEM_OPT_VERSION
   if (opts.defined("binrefcoords"))
   {
     current = config->find("binrefcoords");

     iout << iINFO << "BINARY REF COORDS      "
              << current->data << "\n";
   }
#endif

   if (firstTimestep)
   {
  iout << iINFO << "FIRST TIMESTEP         "
     << firstTimestep << "\n" << endi;
   }
}
/*    END OF FUNCTION initialize_config_data    */


/****************************************************************/
/*                                                              */
/*      FUNCTION parse_mgrid_params                             */
/*                                                              */
/*                                                              */
/****************************************************************/
void SimParameters::parse_mgrid_params(ConfigList *config)
{
  StringList *current;

  mgridforcelist.clear();
  char *key = new char[81];
  char *valstr = new char[256];
  // If the old gridforce commands are still in use, parse them too.
  if (gridforceOn) {
    mgridforceOn = TRUE;
    const char *default_key = MGRIDFORCEPARAMS_DEFAULTKEY;
    MGridforceParams* mgfp = NULL;
    mgfp = mgridforcelist.find_key(default_key);
    if (mgfp != NULL) {
      iout << iINFO << "MGRIDFORCEPOTFILE key "
        << key << " redefined for file " << valstr << "\n" << endi;
    } else {
      mgfp = mgridforcelist.add(default_key);
    }
    mgfp->gridforceVolts = gridforceVolts;
    mgfp->gridforceScale = gridforceScale;

    parse_mgrid_string_param(config,"gridforcefile",&(mgfp->gridforceFile));
    parse_mgrid_string_param(config,"gridforcecol",&(mgfp->gridforceCol));
    parse_mgrid_string_param(config,"gridforcechargecol",&(mgfp->gridforceQcol));
    parse_mgrid_string_param(config,"gridforcepotfile",&(mgfp->gridforceVfile));

    mgfp->gridforceCont[0] = gridforceContA1;
    mgfp->gridforceCont[1] = gridforceContA2;
    mgfp->gridforceCont[2] = gridforceContA3;
    mgfp->gridforceVOffset = gridforceVOffset;

    mgfp->gridforceLite = gridforceLite;
    mgfp->gridforceCheckSize = gridforcechecksize;
  }

  // Create multigrid parameter structures
  current = config->find("mgridforcepotfile");
  while (current != NULL) {
    int curlen = strlen(current->data);
    //    iout << iINFO << "MGRIDFORCEPOTFILE " << current->data
    //         << " " << curlen << "\n"  << endi;
    sscanf(current->data,"%80s%255s",key,valstr);

    MGridforceParams* mgfp = NULL;
    mgfp = mgridforcelist.find_key(key);
    if ( mgfp != NULL) {
      iout << iINFO << "MGRIDFORCEPOTFILE key "
        << key << " redefined for file " << valstr << "\n" << endi;
    } else {
      mgfp = mgridforcelist.add(key);
    }
    int fnamelen = strlen(valstr);
    mgfp->gridforceVfile = new char[fnamelen+1];
    strncpy(mgfp->gridforceVfile,valstr,fnamelen+1);
    mgfp->gridforceScale.x =
      mgfp->gridforceScale.y =
        mgfp->gridforceScale.z = 1.;
    mgfp->gridforceVOffset.x =
      mgfp->gridforceVOffset.y =
        mgfp->gridforceVOffset.z = 0.;

    current = current->next;
  }

  current = config->find("mgridforcefile");
  while (current != NULL) {
    int curlen = strlen(current->data);
    //    iout << iINFO << "MGRIDFORCEFILE " << current->data
    //         << " " << curlen << "\n"  << endi;
    sscanf(current->data,"%80s%255s",key,valstr);

    MGridforceParams* mgfp = NULL;
    mgfp = mgridforcelist.find_key(key);
    if ( mgfp == NULL) {
      iout << iINFO << "MGRIDFORCEFILE no key "
      << key << " defined for file " << valstr << "\n" << endi;
    } else {
      int fnamelen = strlen(valstr);
      if (mgfp->gridforceFile != NULL) {
        delete [] mgfp->gridforceFile;
      }
      mgfp->gridforceFile = new char[fnamelen+1];
      strncpy(mgfp->gridforceFile,valstr,fnamelen+1);
    }

    current = current->next;
  }

  current = config->find("mgridforcevolts");
  while (current != NULL) {
    //    iout << iINFO << "MGRIDFORCEVOLTS " << current->data << "\n"
    //         << endi;
    int curlen = strlen(current->data);
    sscanf(current->data,"%80s%255s",key,valstr);

    MGridforceParams* mgfp = NULL;
    mgfp = mgridforcelist.find_key(key);
    if ( mgfp == NULL) {
      iout << iINFO << "MGRIDFORCEVOLTS no key "
      << key << " defined for file " << valstr << "\n" << endi;
    } else {
      int boolval = MGridforceParamsList::atoBool(valstr);
      if (boolval == -1) {
        iout << iINFO << "MGRIDFORCEVOLTS  key "
          << key << " boolval " << valstr << " badly defined" << endi;
      } else {
        mgfp->gridforceVolts = (boolval == 1);
      }
    }

    current = current->next;
  }

  current = config->find("mgridforcescale");
  while (current != NULL) {
    //    iout << iINFO << "MGRIDFORCESCALE " << current->data
    //         << "\n"  << endi;
    int curlen = strlen(current->data);
    int nread;
    sscanf(current->data,"%80s%n",key,&nread);
    char *val = current->data + nread + 1;

    MGridforceParams* mgfp = NULL;
    mgfp = mgridforcelist.find_key(key);
    if ( mgfp == NULL) {
      iout << iINFO << "MGRIDFORCESCALE no key "
      << key << " defined for vector " << val << "\n" << endi;
    } else {
      mgfp->gridforceScale.set(val);
    }

    current = current->next;
  }

  current = config->find("mgridforcevoff");
  while (current != NULL) {
    //    iout << iINFO << "MGRIDFORCEVOFF " << current->data
    //         << "\n"  << endi;
    int curlen = strlen(current->data);
    int nread;
    sscanf(current->data,"%80s%n",key,&nread);
    char *val = current->data + nread + 1;

    MGridforceParams* mgfp = NULL;
    mgfp = mgridforcelist.find_key(key);
    if ( mgfp == NULL) {
      iout << iINFO << "MGRIDFORCEVOFF no key "
      << key << " defined for vector " << val << "\n" << endi;
    } else {
      mgfp->gridforceVOffset.set(val);
    }

    current = current->next;
  }

  current = config->find("mgridforcecol");
  while (current != NULL) {
    //    iout << iINFO << "MGRIDFORCECOL " << current->data
    //         << "\n"  << endi;
    int curlen = strlen(current->data);
    sscanf(current->data,"%80s%255s",key,valstr);

    MGridforceParams* mgfp = NULL;
    mgfp = mgridforcelist.find_key(key);
    if ( mgfp == NULL) {
      iout << iINFO << "MGRIDFORCECOL no key "
      << key << " defined for file " << valstr << "\n" << endi;
    } else {
      int collen = strlen(valstr);
      if (mgfp->gridforceCol != NULL) {
        delete [] mgfp->gridforceCol;
      }
      mgfp->gridforceCol = new char[collen+1];
      strncpy(mgfp->gridforceCol,valstr,collen+1);
     }

    current = current->next;
  }

  current = config->find("mgridforcechargecol");
  while (current != NULL) {
    //    iout << iINFO << "MGRIDFORCECHARGECOL " << current->data << "\n"
    //         << endi;
    int curlen = strlen(current->data);
    sscanf(current->data,"%80s%255s",key,valstr);

    MGridforceParams* mgfp = NULL;
    mgfp = mgridforcelist.find_key(key);
    if ( mgfp == NULL) {
      iout << iINFO << "MGRIDFORCECHARGECOL no key "
      << key << " defined for file " << valstr << "\n" << endi;
    } else {
      int collen = strlen(valstr);
      if (mgfp->gridforceQcol != NULL) {
        delete [] mgfp->gridforceQcol;
      }
      mgfp->gridforceQcol = new char[collen+1];
      strncpy(mgfp->gridforceQcol,valstr,collen+1);
    }

    current = current->next;
  }

  current = config->find("mgridforcecont1");
  while (current != NULL) {
    //    iout << iINFO << "MGRIDFORCECONT1 " << current->data
    //         << "\n"  << endi;
    int curlen = strlen(current->data);
    sscanf(current->data,"%80s%255s",key,valstr);

    MGridforceParams* mgfp = NULL;
    mgfp = mgridforcelist.find_key(key);
    if ( mgfp == NULL) {
      iout << iINFO << "MGRIDFORCECONT1 no key "
      << key << " defined for file " << valstr << "\n" << endi;
    } else {
      int boolval = MGridforceParamsList::atoBool(valstr);
      if (boolval == -1) {
        iout << iINFO << "MGRIDFORCECONT1  key "
        << key << " boolval " << valstr << " badly defined" << endi;
      } else {
        mgfp->gridforceCont[0] = (boolval == 1);
      }
    }

    current = current->next;
  }

  current = config->find("mgridforcecont2");
  while (current != NULL) {
    //    iout << iINFO << "MGRIDFORCECONT2 " << current->data
    //         << "\n"  << endi;
    int curlen = strlen(current->data);
    sscanf(current->data,"%80s%255s",key,valstr);

    MGridforceParams* mgfp = NULL;
    mgfp = mgridforcelist.find_key(key);
    if ( mgfp == NULL) {
      iout << iINFO << "MGRIDFORCECONT2 no key "
      << key << " defined for file " << valstr << "\n" << endi;
    } else {
      int boolval = MGridforceParamsList::atoBool(valstr);
      if (boolval == -1) {
        iout << iINFO << "MGRIDFORCECONT2  key "
        << key << " boolval " << valstr << " badly defined" << endi;
      } else {
        mgfp->gridforceCont[1] = (boolval == 1);
      }
    }

    current = current->next;
  }
  current = config->find("mgridforcecont3");
  while (current != NULL) {
    //    iout << iINFO << "MGRIDFORCECONT3 " << current->data
    //         << "\n"  << endi;
    int curlen = strlen(current->data);
    sscanf(current->data,"%80s%255s",key,valstr);

    MGridforceParams* mgfp = NULL;
    mgfp = mgridforcelist.find_key(key);
    if ( mgfp == NULL) {
      iout << iINFO << "MGRIDFORCECONT3 no key "
      << key << " defined for file " << valstr << "\n" << endi;
      NAMD_die("MGRIDFORCE error");
    } else {
      int boolval = MGridforceParamsList::atoBool(valstr);
      if (boolval == -1) {
        iout << iINFO << "MGRIDFORCECONT3  key "
        << key << " boolval " << valstr << " badly defined" << endi;
      } else {
        mgfp->gridforceCont[2] = (boolval == 1);
      }
    }

    current = current->next;
  }

  current = config->find("mgridforcelite");
  while (current != NULL) {
    //    iout << iINFO << "MGRIDFORCELITE " << current->data << "\n"
    //         << endi;
    int curlen = strlen(current->data);
    sscanf(current->data,"%80s%255s",key,valstr);

    MGridforceParams* mgfp = NULL;
    mgfp = mgridforcelist.find_key(key);
    if ( mgfp == NULL) {
      iout << iINFO << "MGRIDFORCELITE no key "
      << key << " defined for file " << valstr << "\n" << endi;
    } else {
      int boolval = MGridforceParamsList::atoBool(valstr);
      if (boolval == -1) {
        iout << iINFO << "MGRIDFORCELITE  key "
          << key << " boolval " << valstr << " badly defined" << endi;
      } else {
        mgfp->gridforceLite = (boolval == 1);
      }
    }

    current = current->next;
  }

  current = config->find("mgridforcechecksize");
  while (current != NULL) {
    //    iout << iINFO << "MGRIDFORCELITE " << current->data << "\n"
    //         << endi;
    int curlen = strlen(current->data);
    sscanf(current->data,"%80s%255s",key,valstr);

    MGridforceParams* mgfp = NULL;
    mgfp = mgridforcelist.find_key(key);
    if ( mgfp == NULL) {
      iout << iINFO << "MGRIDFORCECHECKSIZE no key "
      << key << " defined for file " << valstr << "\n" << endi;
    } else {
      int boolval = MGridforceParamsList::atoBool(valstr);
      if (boolval == -1) {
        iout << iINFO << "MGRIDFORCECHECKSIZE  key "
          << key << " boolval " << valstr << " badly defined" << endi;
      } else {
        mgfp->gridforceCheckSize = (boolval == 1);
      }
    }

    current = current->next;
  }

  delete [] valstr;
  delete [] key;

  // Fill in default values for optional items

  MGridforceParams* params = mgridforcelist.get_first();

  while (params != NULL) {
    if (params->gridforceFile == NULL) {
      char errmsg[255];
      sprintf(errmsg,"Value undefined for gridforceFile for key %s\n",
              params->gridforceKey);
         NAMD_die(errmsg);
    }
    if (params->gridforceCol == NULL) {
      char errmsg[255];
      sprintf(errmsg,"Value undefined for gridforceCol for key %s\n",
              params->gridforceKey);
         NAMD_die(errmsg);
    }
    params = params->next;
  }

}

void SimParameters::parse_mgrid_string_param(ConfigList *cl,
                                             const char *fieldname,
                                             char **dest)
{
  StringList *vallist = cl->find(fieldname);
  char *val = NULL;

  if (vallist != NULL) {
    val = vallist->data;
  } else {
    return;
  }

  int len = 0;
  if (val == NULL) {
    *dest = NULL;
  } else {
    len = strlen(val);
    if (len == 0) {
      *dest = NULL;
    } else {
      *dest = new char[len+1];
      strncpy(*dest,val,len+1);
    }
  }
}

//This function is used to create directories when outputing into
//multiple files, i.e. used for Parallel IO. -Chao Mei
void SimParameters::create_output_directories(const char *dirname){
	//output files organization:
	//$outputFilename/$dirname/$outputproc_rank

	//Step 1: create $outputFilename if necessary
	int baselen = strlen(outputFilename);
	char *filename = new char[baselen+32];
	memset(filename, 0, baselen+32);
	strcpy(filename, outputFilename);
	if(access(filename, F_OK)!=0) {
		int ret = MKDIR(filename);
		if(ret!=0) {
			char errmsg[512];
			sprintf(errmsg, "Error in creating top-level directory %s!", filename);
			NAMD_die(errmsg);
		}
	}

	//Step 2: create $dirname if necessary
	strcat(filename, PATHSEPSTR);
	strcat(filename, dirname);
	//check if the directory exists or not
	if(access(filename, F_OK)!=0) {
		int ret = MKDIR(filename);
		if(ret!=0) {
			char errmsg[512];
			sprintf(errmsg, "Error in creating middle-level directory %s!", filename);
			NAMD_die(errmsg);
		}
	}

	//step 3: create $outputproc_rank if necessary
	char tmpstr[256];
	for(int i=0; i<numoutputprocs; i++) {
		memset(tmpstr, 0, 256);
		sprintf(tmpstr, "%s%s%d", filename, PATHSEPSTR, i);
		if(access(tmpstr, F_OK)!=0) {
			int ret = MKDIR(tmpstr);
			if(ret!=0) {
				char errmsg[512];
				sprintf(errmsg, "Error in creating last-level directory %s!", tmpstr);
				NAMD_die(errmsg);
			}
		}
	}
}

/****************************************************************/
/*                                                              */
/*      FUNCTION print_mgrid_params                             */
/*                                                              */
/*                                                              */
/****************************************************************/
#define BoolToString(b) ((b) ? "TRUE" : "FALSE")

void SimParameters::print_mgrid_params()
{
  const MGridforceParams* params = mgridforcelist.get_first();

  while (params != NULL) {
    iout << iINFO << "MGRIDFORCE key " << params->gridforceKey << "\n" << endi;
    iout << iINFO << "           Potfile " << params->gridforceVfile
      << "\n" << endi;
    iout << iINFO << "           Scale " << params->gridforceScale
      << "\n" << endi;
    iout << iINFO << "           File " << params->gridforceFile
      << "\n" << endi;
    iout << iINFO << "           Col " << params->gridforceCol
      << "\n" << endi;

    const char *qcol_msg = "Use atom charge";
    if (params->gridforceQcol != NULL) {
      qcol_msg = params->gridforceQcol;
    }
    iout << iINFO << "           ChargeCol " << qcol_msg
      << "\n" << endi;
    iout << iINFO << "           VOffset " << params->gridforceVOffset
      << "\n" << endi;
    iout << iINFO << "           Continuous K1 "
      << BoolToString(params->gridforceCont[0])
      << "\n" << endi;
    iout << iINFO << "           Continuous K2 "
      << BoolToString(params->gridforceCont[1])
    << "\n" << endi;
    iout << iINFO << "           Continuous K3 "
      << BoolToString(params->gridforceCont[2])
      << "\n" << endi;
    iout << iINFO << "           Volts "
      << BoolToString(params->gridforceVolts)
      << "\n" << endi;
    iout << iINFO << "           Gridforce-Lite "
      << BoolToString(params->gridforceLite)
      << "\n" << endi;
    iout << iINFO << "           Gridforce-CheckSize "
      << BoolToString(params->gridforceCheckSize)
      << "\n" << endi;
    params = params->next;
  }
}


/****************************************************************/
/*                */
/*    FUNCTION send_SimParameters      */
/*                */
/*  This function is used by the master process to broadcast*/
/*  the parameter data to all the other nodes.  It just builds  */
/*  a message with all the relevant data and broadcasts it to   */
/*  the other nodes.  The routine receive_SimParameters is used */
/*  by all the other nodes to receive this message.    */
/*                */
/****************************************************************/

void SimParameters::send_SimParameters(MOStream *msg)

{
  /*MOStream *msg = com_obj->newOutputStream(ALLBUTME, SIMPARAMSTAG, BUFSIZE);
  if ( msg == NULL )
  {
    NAMD_die("memory allocation failed in SimParameters::send_SimParameters");
  }*/

  msg->put(sizeof(SimParameters),(char*)this);
  if ( FFTWWisdomString ) {
    int fftwlen = strlen(FFTWWisdomString) + 1;
    msg->put(fftwlen);
    msg->put(fftwlen,FFTWWisdomString);
  }
  if ( tclBCScript ) {
    int tcllen = strlen(tclBCScript) + 1;
    msg->put(tcllen);
    msg->put(tcllen,tclBCScript);
  }

#ifdef MEM_OPT_VERSION
  int filelen = strlen(binAtomFile)+1;
  msg->put(filelen);
  msg->put(filelen, binAtomFile);

  filelen = strlen(binCoorFile)+1;
  msg->put(filelen);
  msg->put(filelen, binCoorFile);

  if(binVelFile) {
    filelen = strlen(binVelFile)+1;
    msg->put(filelen);
    msg->put(filelen, binVelFile);
  }

  if(binRefFile) {
    filelen = strlen(binRefFile)+1;
    msg->put(filelen);
    msg->put(filelen, binRefFile);
  }
#endif

  mgridforcelist.pack_data(msg);

  msg->end();
}
/*    END OF FUNCITON send_SimParameters    */

/****************************************************************/
/*                */
/*      FUNCTION receive_SimParameters    */
/*                */
/*  This function is used by all the child nodes to   */
/*  receive the simulation parameters from the master node.  */
/*                */
/****************************************************************/

void SimParameters::receive_SimParameters(MIStream *msg)

{
  msg->get(sizeof(SimParameters),(char*)this);
  if ( FFTWWisdomString ) {
    int fftwlen;
    msg->get(fftwlen);
    FFTWWisdomString = new char[fftwlen];
    msg->get(fftwlen,FFTWWisdomString);
#ifdef NAMD_FFTW
#ifdef NAMD_FFTW_3
    fftwf_import_wisdom_from_string(FFTWWisdomString);
#else
    fftw_import_wisdom_from_string(FFTWWisdomString);
#endif
#endif
  }
  if ( tclBCScript ) {
    int tcllen;
    msg->get(tcllen);
    tclBCScript = new char[tcllen];
    msg->get(tcllen,tclBCScript);
  }

#ifdef MEM_OPT_VERSION
  int filelen;
  msg->get(filelen);
  binAtomFile = new char[filelen];
  msg->get(filelen, binAtomFile);

  msg->get(filelen);
  binCoorFile = new char[filelen];
  msg->get(filelen, binCoorFile);

  if(binVelFile) {
    msg->get(filelen);
    binVelFile = new char[filelen];
    msg->get(filelen, binVelFile);
  }

  if(binRefFile) {
    msg->get(filelen);
    binRefFile = new char[filelen];
    msg->get(filelen, binRefFile);
  }
#endif


  // The simParameters bit copy above put illegal values in the list pointers
  // So this resets everything so that unpacking will work.
  mgridforcelist.clear();
  mgridforcelist.unpack_data(msg);

  delete msg;
}
/*      END OF FUNCTION receive_SimParameters  */


//fepb IDWS
BigReal SimParameters::getCurrentLambda2(const int step) {
  if ( alchLambdaIDWS >= 0. ) {
    const BigReal lambda2 = ( (step / alchIDWSFreq) % 2 == 1 ) ? alchLambda2 : alchLambdaIDWS;
    return lambda2;
  } else {
    return alchLambda2;
  }
}

/* Return true if IDWS is active, else return false. */
int SimParameters::setupIDWS() {
  if (alchLambdaIDWS < 0.) return 0;
  if (alchLambdaIDWS > 1.) {
    NAMD_die("alchLambdaIDWS should be either in the range [0.0, 1.0], or negative (disabled).\n");
  }
 /*
  * The internal parameter alchIDWSFreq determines the number of steps of MD
  * before each switch of the value of alchLambda2. At most this occurs every
  * time the energy is evaluated and thus the default is the greater of
  * fullElectFrequency and nonbondedFrequency. However, this choice fails to
  * report alternating values if output is printed less often than every step
  * (which is almost certainly true). Thus the frequency is reset to match
  * alchOutFreq or, if that is zero, outputEnergies. Note that, if
  * alchOutFreq > 0 but != outputEnergies, then the data going to stdout
  * are likely not useful since the comparison value is difficult to infer.
  */
  alchIDWSFreq = fullElectFrequency > 0 ? fullElectFrequency : nonbondedFrequency;
  if ( !alchOutFreq && outputEnergies > alchIDWSFreq ) alchIDWSFreq = outputEnergies;
  if ( alchOutFreq > alchIDWSFreq ) alchIDWSFreq = alchOutFreq;
  if ( alchOutFreq && alchOutFreq != outputEnergies) {
    iout << iWARN << "alchOutFreq and outputEnergies do not match. IDWS ouput"
         << " to stdout may not be useful!\n" << endi;
  }
  return 1;
}
//fepe IDWS

//fepb BKR
BigReal SimParameters::getCurrentLambda(const int step) {
  /*Get lambda at the current step.

   If alchLambdaFreq = 0, return alchLambda. For positive values of
   alchLambdaFreq, apply a linear stepwise schedule from alchLambda to
   alchLambda2:

   l(t) = l + (l2 - l)*[dn / (N - n0)]*{floor[(n - n0)/dn] + 1}

   n - the current time step
   n0 - step at which switching begins (default = 0)
   N - total steps in the simulation
   dn - alchLambdaFreq (increment frequency, in steps)
   l/l2 - alchLambda/alchLambda2

   Note that each step _begins_ by incrementing alchLambda and then integrates
   in time. This means that the first and last switch steps may not behave as
   immediately expected - at step 0, alchLambda is NOT evaluated and at step N
   no step occurs because alchLambda2 has already been reached.
  */
  if ( alchLambdaFreq > 0 && step >= alchEquilSteps ) {
    if ( step == N ) {
      return alchLambda2;
    }
    else {
      const int timeOrigin = firstTimestep + alchEquilSteps;
      const BigReal alchLambdaDelta = getLambdaDelta();
      const BigReal increment = (step - timeOrigin) / BigReal(alchLambdaFreq);
      return alchLambda + alchLambdaDelta*(floor(increment) + 1);
    }
  }
  else {
    return alchLambda;
  }
}

BigReal SimParameters::getLambdaDelta(void) {
  // Increment by which Lambda changes.
  return ((alchLambda2 - alchLambda)*alchLambdaFreq
          / BigReal(N - firstTimestep - alchEquilSteps));
}

BigReal SimParameters::getElecLambda(const BigReal lambda) {
  // Convenience function for staggered lambda scaling
  return (lambda <= alchElecLambdaStart ? 0.
          : (lambda - alchElecLambdaStart) / (1. - alchElecLambdaStart));
}

/*
 * Modifications for WCA decomposition of van der Waal interactions.
 *
 * WCA requires that repulsive and attractive components of the vdW
 * forces be treated separately. To keep the code clean, the same scaling
 * function is always used and simply has its behavior modified. However,
 * the new repluslive scaling only ever gets used when alchWCAOn.
 */
BigReal SimParameters::getVdwLambda(const BigReal lambda) {
  // Convenience function for staggered lambda scaling
  if ( alchWCAOn ) {
    // Read this with the alias alchRepLambdaEnd --> alchAttLambdaStart.
    // The second condition is needed when attractive interactions are inactive
    // for the whole range, otherwise lambda = 0/1 are incorrect.
    if ( lambda < alchRepLambdaEnd || alchRepLambdaEnd == 1.0 ) {
      return 0.0;
    } else if ( lambda >= alchVdwLambdaEnd ) {
      return 1.0;
    } else {
      return (lambda - alchRepLambdaEnd) / (alchVdwLambdaEnd - alchRepLambdaEnd);
    }
  } else {
  return (lambda >= alchVdwLambdaEnd ? 1. : lambda / alchVdwLambdaEnd);
}
}

BigReal SimParameters::getRepLambda(const BigReal lambda) {
  // Convenience function for staggered lambda scaling
  return (lambda >= alchRepLambdaEnd ? 1. : lambda / alchRepLambdaEnd);
}

BigReal SimParameters::getBondLambda(const BigReal lambda) {
  // Convenience function for staggered lambda scaling
  return (lambda >= alchBondLambdaEnd ? 1. : lambda / alchBondLambdaEnd);
}
//fepe
