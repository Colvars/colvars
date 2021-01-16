/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*****************************************************************************
 * $Source: /home/cvs/namd/cvsroot/namd2/src/SimParameters.h,v $
 * $Author: jim $
 * $Date: 2017/03/30 20:06:17 $
 * $Revision: 1.1248 $
 *****************************************************************************/

#ifndef SIMPARAMETERS_H
#define SIMPARAMETERS_H

#include "common.h"
#include "Vector.h"
#include "Lattice.h"

#include "MGridforceParams.h"

class ParseOptions;
class Communicate;
class ConfigList;
class MIStream;

//  The class SimParameters is really just a glorified structure used to
//  maintain the global simulation parameters.  The only functions
//  associated with the class are used to get the parameters from the
//  ConfigList object, to send that Parameters from the master node
//  to the other nodes, and to receive the Parameters on the other nodes.


//  The following definitions are used to distinguish between possible
//  bonded exclusion settings
typedef int  ExclusionSettings;

#define NONE		0
#define	ONETWO  	1
#define	ONETHREE	2
#define	ONEFOUR		3
#define	SCALED14 	4

//  The following definitions are used to distinguish between multiple
//  timestep integration schemes
typedef int  MTSChoices;

#define NAIVE		0
#define VERLETI		1

//  The following definitions are used to distinuish between multiple
//  long-short range force splittings
#define SHARP		0
#define XPLOR		1
#define C1		2
#define C2              3

//  The following definitions are used to distinguish among load
//  balancers and their strategies
#define LDBAL_NONE		0
#define LDBAL_CENTRALIZED	1	// default
#define LDBAL_HYBRID		2

#define LDBSTRAT_DEFAULT	10	// default
#define LDBSTRAT_COMPREHENSIVE	11
#define LDBSTRAT_REFINEONLY	12
#define LDBSTRAT_OLD		13

// The following definitions are used to distinguish between patch-splitting
// strategies
#define SPLIT_PATCH_POSITION	0	// atom position determines patch
#define SPLIT_PATCH_HYDROGEN	1	// hydrogen groups are not broken up

// The following definitions are used to distinguish the range of rigid
// bond calculations: none, all bonds to hydrogen, or only water
#define RIGID_NONE    0
#define RIGID_ALL     1
#define RIGID_WATER   2

// Added by JLai -- The following definitions are used to distinguish
// the different GoMethodologies available to the Go program
// -- 6.3.11
typedef int GoChoices;
#define GO_MATRIX 1
#define GO_SPARSE 2
#define GO_LOWMEM 3

// Used for controlling PME parallelization with ckloop
// The higher level will include all parallelization for lower ones
// E.g. If setting useCkLoop to 3, then xpencil's kspace, all
// backward ffts and send_untrans/ungrid routines will be parallelized
#define CKLOOP_CTRL_PME_UNGRIDCALC 6
#define CKLOOP_CTRL_PME_FORWARDFFT 5
#define CKLOOP_CTRL_PME_SENDTRANS 4
#define CKLOOP_CTRL_PME_KSPACE 3
#define CKLOOP_CTRL_PME_BACKWARDFFT 2
#define CKLOOP_CTRL_PME_SENDUNTRANS 1

class SimParameters
{
private:
public:

//  MAKE SURE THAT THIS CLASS CAN BE BIT COPIED OR YOU WILL HAVE TO
//  ADD SPECIAL CODE TO send_SimParameters() and receive_SimParameters()

#if defined(NAMD_NVTX_ENABLED) || defined(NAMD_CMK_TRACE_ENABLED)
  int beginEventPatchID;
  int endEventPatchID;
  int beginEventStep;
  int endEventStep;
#endif

#ifdef TIMER_COLLECTION
  double timerBinWidth;  // default 1
#endif

  Bool lonepairs;  // enable lone pairs
  int watmodel; // integer code for the water model in use
                // choices are defined in common.h
  Bool LJcorrection; // flag for whether water tail corrections should be used
	BigReal dt;	   		//  Timestep size
	int N;		   		//  Number of steps to be performed
	int stepsPerCycle;		//  Number of timesteps per cycle

	zVector cellBasisVector1;	//  Basis vector for periodic cell
	zVector cellBasisVector2;	//  Basis vector for periodic cell
	zVector cellBasisVector3;	//  Basis vector for periodic cell
	zVector cellOrigin;		//  Fixed center of periodic cell
	Lattice lattice;		//  All data for periodic cell

	int nonbondedFrequency;		//  Number of timesteps between
					//  nonbonded evaluation
	int fullElectFrequency;		//  Number of timesteps between
					//  full electrostatic evaluation
        BigReal fmaTheta;	        //  DPMTA theta value
	int ldBalancer;			//  None, Centralized or Hybrid
	int ldbStrategy;                //  What load balancing strategy to use
	int ldbPeriod;                  //  How often to do load balancing
	int firstLdbStep;		//  What step to do the first
                                        //  load-balance on.
	int lastLdbStep;		//  What step to do the last
                                        //  load-balance on.
        int hybridGroupSize;            //  hybrid group size
	BigReal ldbBackgroundScaling;	//  scaling factor for background load
	BigReal ldbPMEBackgroundScaling;//  scaling factor for PME background
	BigReal ldbHomeBackgroundScaling;//  scaling factor for home background
	BigReal ldbRelativeGrainsize;   //  fraction of average load per compute

	int traceStartStep; //the timestep when trace is turned on, default to 3*firstLdbStep;
	int numTraceSteps; //the number of timesteps that are traced, default to 2*ldbPeriod;

#ifdef MEASURE_NAMD_WITH_PAPI
	Bool papiMeasure; //default to false
	int papiMeasureStartStep; //the timestep when to measure using PAPI, default to 3*firstLdbStep;
	int numPapiMeasureSteps; //the number of timesteps when performance are measured with PAPI, default to 40;
#endif

	Bool outputMaps; //control whether to dump compute/patch map before load balancing
	Bool simulateInitialMapping; //if true, the initial mapping during startup is dumped and exit
	int simulatedPEs;
	int simulatedNodeSize;
	Bool disableTopology; // ignore torus information during patch placement
	Bool verboseTopology; // print torus information during patch placement

	Bool benchTimestep; //only cares about benchmarking the timestep, so no file output to save SUs for large-scale benchmarking

	//whether to use CkLoop library to parallelize a loop in a function like OpenMP.
	//It has multiple control levels. The higher the value is (must be positive), the more parallelization will be performed
	//Currently, it is mainly used for PME computation. The default value is 0, meaning it is disabled
	//Refer to macros CKLOOP_CTRL_* in this file for the ordering of different levels
	int useCkLoop;

	int twoAwayX;			//  half-size patches in X dimension
	int twoAwayY;			//  half-size patches in Y dimension
	int twoAwayZ;			//  half-size patches in Z dimension
	int maxPatches;			//  maximum patch count
	Bool ldbUnloadPME;		//  unload processors doing PME
	Bool ldbUnloadZero;		//  unload processor 0
	Bool ldbUnloadOne;		//  unload processor 1
	Bool ldbUnloadOutputPEs;	//  unload output processors
	Bool noPatchesOnZero;		//  no patches on processor 0
	Bool noPatchesOnOutputPEs;	//  no patches on output PEs
	Bool noPatchesOnOne;		//  no patches on processor 1

	BigReal initialTemp;   		//  Initial temperature for the
					//  simulation
	Bool comMove;     		//  Should the center of mass be
					//  able to move
	Bool zeroMomentum;		//  remove momentum drift from PME
	Bool zeroMomentumAlt;		//  alternate method for testing
	Bool wrapWater;			//  Wrap water around on output
	Bool wrapAll;			//  Wrap clusters around on output
	Bool wrapNearest;		//  Wrap to closest image to origin
	BigReal dielectric;   		//  Dielectric constant
	ExclusionSettings exclude;      //  What electrostatic exclusions should
					//  be made
	BigReal scale14;		//  Scaling factor for 1-4
					//  electrostatics
	BigReal nonbondedScaling;	//  Scaling factor for nonbonded forces
	int dcdFrequency;		//  How often (in timesteps) should
					//  a DCD trajectory file be updated
  int dcdUnitCell;  // Whether to write unit cell information in the DCD
	int velDcdFrequency;		//  How often (in timesteps) should
					//  a velocity DCD file be updated
	int forceDcdFrequency;		//  How often (in timesteps) should
					//  a force DCD file be updated
	int xstFrequency;		//  How often (in timesteps) should
					//  a XST trajectory file be updated
	char auxFilename[NAMD_FILENAME_BUFFER_SIZE];		//  auxilary output filename
	char dcdFilename[NAMD_FILENAME_BUFFER_SIZE];		//  DCD filename
	char velDcdFilename[NAMD_FILENAME_BUFFER_SIZE];       //  Velocity DCD filename
	char forceDcdFilename[NAMD_FILENAME_BUFFER_SIZE];     //  Force DCD filename
	char xstFilename[NAMD_FILENAME_BUFFER_SIZE];		//  Extended system trajectory filename
	char outputFilename[NAMD_FILENAME_BUFFER_SIZE];	//  Output file name.  This name will
					//  have .coor appended to it
					//  for the coordinates and
					//  .vel appended to
					//  it for the velocities
	char restartFilename[NAMD_FILENAME_BUFFER_SIZE];	//  Base name of the restart file
	int restartFrequency;		//  How often (in timesteps) shoud the
					//  restart files be updated
        Bool restartSave;		//  unique filenames for restart files
        Bool restartSaveDcd;		//  unique filenames for DCD files
	Bool binaryRestart;		//  should restart files be
					//  binary format rather than PDB
	Bool binaryOutput;		//  should output files be
					//  binary format rather than PDB
	BigReal cutoff;			//  Cutoff distance
	BigReal margin;			//  Fudge factor on patch size
	BigReal patchDimension;		//  Dimension of each side of a patch
					//  This is either cutoff+margin or
					//  pairlistDist+margin depending on
					//  whether or not switching is on
					//  or not
	BigReal limitDist;		//  Distance below which nonbonded
					//  forces between atoms are limited
	Bool switchingActive;		//  Flag TRUE->using switching function
					//  for electrostatics and vdw
	Bool vdwForceSwitching;		//  Flag TRUE->using force switching
					//  function for vdw
	BigReal switchingDist;		//  Distance at which switching
					//  becomes active
	Bool martiniSwitching;		//  Flag TRUE->use Martini residue-based
                                        //  coarse-grain switching function
	Bool martiniDielAllow;          //  Allow non-standard dielectric constant
                                        //  for use with Martini when dielectric != 15.0
	BigReal pairlistDist;		//  Distance within which atom pairs
					//  should be added to pairlist
	int pairlistMinProcs;		//  Minimum number of processors
					//  to enable pairlists
	int usePairlists;		//  Derived from pairlistMinProcs

	int pairlistsPerCycle;		//  regenerate x times per cycle
	BigReal pairlistShrink;		//  tol *= (1 - x) on regeneration
	BigReal pairlistGrow;		//  tol *= (1 + x) on trigger
	BigReal pairlistTrigger;	//  trigger is atom > (1 - x) * tol
	int outputPairlists;		//  print pairlist warnings this often

	Bool constraintsOn;		//  Flag TRUE-> harmonic constraints
					//  active
	int constraintExp;		//  Exponent for harmonic constraints

	/* BEGIN gf */
	Bool gridforceOn;		//  Flag TRUE -> gridforce active
	Bool gridforceVolts;		//  Flag TRUE -> gridforce using volts as units
	zVector gridforceScale;		//  Gridforce scale factor
	Bool gridforceContA1;		//  Flag TRUE -> grid continuous in A1 direction
	Bool gridforceContA2;		//  Flag TRUE -> grid continuous in A2 direction
	Bool gridforceContA3;		//  Flag TRUE -> grid continuous in A3 direction
	zVector gridforceVOffset;	//  Gridforce potential offsets
	Bool gridforceLite;		//  Flag TRUE -> use lightweight, fast, feature-poor gridforce
	Bool gridforcechecksize; //Flag TRUE -> check if grid is larger than PBC cell dimensions
  /* END gf */
        Bool mgridforceOn;
        MGridforceParamsList mgridforcelist;

        //****** BEGIN selective restraints (X,Y,Z) changes
        Bool selectConstraintsOn;       //  Flag TRUE-> selective restraints
                                        //  active
        Bool constrXOn, constrYOn,
             constrZOn;                 //  Flag TRUE-> select which Cartesian
                                        //  component to restrain
        //****** END selective restraints (X,Y,Z) changes

	// spherical constraints
	Bool sphericalConstraintsOn;
	zVector sphericalConstrCenter;

	BigReal constraintScaling;	//  Scaling factor for constraint forces

        //****** BEGIN CHARMM/XPLOR type changes
        Bool paraTypeXplorOn;           //  FLAG TRUE-> parametrs are XPLOR format (default)
        Bool paraTypeCharmmOn;          //  FLAG TRUE-> parametrs are CHARMM format
        //****** END CHARMM/XPLOR type changes

	// Ported by JLai -- JE - Go
        Bool goGroPair;           //  FLAG FALSE->Explicit Gromacs pairs will be calculated
        Bool goForcesOn;          //  FLAG TRUE-> Go forces will be calculated
        char goParameters[NAMD_FILENAME_BUFFER_SIZE];   //  File for Go parameters
        char goCoordinates[NAMD_FILENAME_BUFFER_SIZE];  //  File for Go structure and atom chain types
        //JLai 6.3.11
	GoChoices  goMethod;      //  Integer for Go method -- 1) Matrix-Go, 3) Low-mem-Go
	// End of port -- JL

        //****** BEGIN moving constraints changes
        Bool movingConstraintsOn;       //  Flag TRUE-> moving constraints
                                        //  active
        zVector movingConsVel;           //  Velocity of the movement, A/timestep
        //****** END moving constraints changes
        //****** BEGIN rotating constraints changes
        Bool rotConstraintsOn;          //  Flag TRUE-> rotating constraints
                                        //  active
        zVector rotConsAxis;             //  Axis of rotation
        zVector rotConsPivot;            //  Pivot point of rotation
        BigReal rotConsVel;             //  Velocity of rotation, Deg/timestep
        //****** END rotating constraints changes

        //****** BEGIN moving drag changes
        Bool movDragOn;               //  Flag TRUE-> moving drag active
        char movDragFile[NAMD_FILENAME_BUFFER_SIZE];        //  PDB file defining dragged atoms
                                      //  by non-zero value in the column
	BigReal movDragGlobVel;       //  global drag velocity (A/step)
	char movDragVelFile[NAMD_FILENAME_BUFFER_SIZE];     //  PDB file; XYZ scale moving drag
                                      //  velocity for each atom
        //****** END moving drag changes
        //****** BEGIN rotating drag changes
        Bool rotDragOn;               //  Flag TRUE-> rotating drag active
        char rotDragFile[NAMD_FILENAME_BUFFER_SIZE];        //  PDB file defining dragged atoms
                                      //  by non-zero value in the column
	char rotDragAxisFile[NAMD_FILENAME_BUFFER_SIZE];    //  PDB file; XYZ define axes for atoms;
	char rotDragPivotFile[NAMD_FILENAME_BUFFER_SIZE];   //  PDB file; XYZ define pivots for atoms
	BigReal rotDragGlobVel;       //  global drag velocity (deg/step)
	char rotDragVelFile[NAMD_FILENAME_BUFFER_SIZE];     //  PDB file; B or O scales angular
                                      //  velocity for each atom
        //****** END rotating drag changes
        //****** BEGIN "constant" torque changes
        Bool consTorqueOn;            //  Flag TRUE-> "constant" torque active
        char consTorqueFile[NAMD_FILENAME_BUFFER_SIZE];     //  PDB file defining torqued atoms
                                      //  by non-zero value in the column
	char consTorqueAxisFile[NAMD_FILENAME_BUFFER_SIZE]; //  PDB file; XYZ define axes for atoms;
	char consTorquePivotFile[NAMD_FILENAME_BUFFER_SIZE];//  PDB file; XYZ define pivots for atoms
	BigReal consTorqueGlobVal;    //  global "torque" (Kcal/(mol*A^2))
	char consTorqueValFile[NAMD_FILENAME_BUFFER_SIZE];  //  PDB file; B or O scales "torque"
                                      //  for each atom
        //****** END "constant" torque changes

        //****** BEGIN SMD constraints changes
        Bool SMDOn;                     //  Flag TRUE-> SMD constraints active
        BigReal SMDVel;                 //  Velocity of the movement, A/timestep
        zVector SMDDir;                  //  Direction of the movement
        BigReal SMDk; 			//  Elastic constant for SMD
	BigReal SMDk2;			//  Transverse elastic constant for SMD
 	char SMDFile[NAMD_FILENAME_BUFFER_SIZE];		//  File for SMD information
        int SMDOutputFreq;              //  Output frequency for SMD constr.
        //****** END SMD constraints changes

  //****** BEGIN tabulated energy section
  Bool tabulatedEnergies;
  int tableNumTypes;
  char tabulatedEnergiesFile[NAMD_FILENAME_BUFFER_SIZE];
  char tableInterpType[128];
  Real tableSpacing;
  BigReal tableMaxDist;
  //****** END tabulated energy section

        // TMD
        Bool TMDOn, TMDDiffRMSD;
        BigReal TMDk;
        char TMDFile[NAMD_FILENAME_BUFFER_SIZE], TMDFile2[NAMD_FILENAME_BUFFER_SIZE];
        int TMDOutputFreq;
        int TMDFirstStep, TMDLastStep;
        BigReal TMDInitialRMSD, TMDFinalRMSD;

        //Symmetry restraints
        Bool symmetryOn, symmetryScaleForces;
        BigReal symmetryk;
        char symmetrykfile[NAMD_FILENAME_BUFFER_SIZE];
        char symmetryFile[NAMD_FILENAME_BUFFER_SIZE];
        char symmetryMatrixFile[NAMD_FILENAME_BUFFER_SIZE];
        int symmetryFirstStep, symmetryLastStep, symmetryFirstFullStep, symmetryLastFullStep;


//fepb
  Bool alchOnAtStartup;     //  Ensure that alchemy is set up properly
  Bool alchFepOnAtStartup;
  Bool alchThermIntOnAtStartup;
  Bool alchOn;              //  Doing alchemical simulation?
  Bool alchFepOn;           //  Doing alchemical simulation?
  Bool singleTopology;      //  Using single topology setup? 
  Bool sdScaling;           //  Scaling S-D bond terms in single topology?
  Bool alchThermIntOn;      //  Doing thermodynamic integration?
  Bool alchWCAOn;           //  Using WCA decomposition for vdWs?
  int alchMethod;           //  Which alchemical method to use? fep or ti
  BigReal alchLambda;       //  lambda for dynamics
  BigReal alchLambda2;      //  lambda for comparison
  BigReal alchLambdaIDWS;   //  alternate lambda for interleaved double-wide sampling
  int alchIDWSFreq;         //  freq with which lambda2 changes to lambdaIDWS
  int alchLambdaFreq;       //  freq. (in steps) with which lambda changes
                            //  from alchLambda to alchLambda2
  BigReal getCurrentLambda(const int); // getter for changing lambda
  BigReal getCurrentLambda2(const int); // getter for alternating lambda2 in IDWS
  int setupIDWS();          //  activates IDWS and sets alchIDWSFreq
  BigReal getLambdaDelta(void); // getter for lambda increment
  BigReal alchTemp;         //  temperature for alchemical calculation
  int alchOutFreq;          //  freq. of alchemical output
  Bool alchEnsembleAvg;      //if do ensemble average for the net free energy difference
  char alchOutFile[NAMD_FILENAME_BUFFER_SIZE];    //  alchemical output filename
  int alchEquilSteps;       //  # of equil. steps in the window
  BigReal alchVdwShiftCoeff; //  r2 shift coeff used for generating
                            //  the alchemical altered vdW interactions
  BigReal alchElecLambdaStart;  //  lambda value for starting point of
                                //  electrostatic interactions of
                                //  exnihilated particles.  For annihilated
                                //  particles the starting point is
                                //  (1-alchElecLambdaStart)
  BigReal getElecLambda(const BigReal); // return min[0,x/(1-elecStart)]
  BigReal alchVdwLambdaEnd;  //  lambda value for endpoint of vdW
                             //  interactions of exnihilated particles.
                             //  For annihilated particles the endpoint is
                             //  (1-alchVdwLambdaEnd)
  BigReal getVdwLambda(const BigReal); // return max[1,x/vdwEnd]
  BigReal alchRepLambdaEnd;  //  lambda value for endpoint of repulsive vdW
                             //  interactions of exnihilated particles.
                             //  For annihilated particles the endpoint is
                             //  (1-alchRepLambdaEnd). This also implies the
                             //  START for attractive vdW interactions.
  BigReal getRepLambda(const BigReal); // return max[1,x/repEnd]
  BigReal alchBondLambdaEnd; //  lambda value for endpoint of bonded
                             //  interactions involving exnihilated particles.
                             //  For annihilated particles the endpoint is
                             //  (1-alchBondLambdaEnd)
  BigReal getBondLambda(const BigReal); // return max[1,x/bondEnd]
  Bool alchDecouple;  // alchemical decoupling rather than annihilation
  Bool alchBondDecouple; // decouple purely alchemical bonds
//fepe


	Bool lesOn;			//  Locally enhanced sampling?
	int lesFactor;			//  local enhancement factor
	Bool lesReduceTemp;		//  Reduce enhanced atom temperature?
	Bool lesReduceMass;		//  Reduce enhanced atom mass?

  // REST2
  Bool soluteScalingOn;
  /**< REST2 (replica exchange solute tempering) on? */
  BigReal soluteScalingFactor;
  /**< scaling factor for solute interactions */
  BigReal soluteScalingFactorCharge;
  /**< optional independent control over scaling factor for
   * electrostatic interactions of solute */
  BigReal soluteScalingFactorVdw;
  /**< optional independent control over scaling factor for
   * van der Waals interactions of solute */
  Bool soluteScalingAll;
  /**< enables scaling for bond and angle terms (default is off) */

        Bool extForcesOn;		//  Are ext command forces present?
        char extForcesCommand[NAMD_FILENAME_BUFFER_SIZE];
        char extCoordFilename[NAMD_FILENAME_BUFFER_SIZE];
        char extForceFilename[NAMD_FILENAME_BUFFER_SIZE];


        // Defines variables for QM/MM calculations
        Bool qmForcesOn;               //  Are QM/MM command forces present?
        char qmParamPDB[NAMD_FILENAME_BUFFER_SIZE];
        Bool qmParamPDBDefined;
        Bool qmChrgFromPSF;
        char qmExecPath[NAMD_FILENAME_BUFFER_SIZE];
        char qmSoftware[128];
        char qmChrgModeS[16];
        int qmChrgMode;
        char qmColumn[16];
        char qmBaseDir[NAMD_FILENAME_BUFFER_SIZE];
        char qmSecProc[NAMD_FILENAME_BUFFER_SIZE];
        Bool qmSecProcOn;
        char qmPrepProc[NAMD_FILENAME_BUFFER_SIZE];
        Bool qmPrepProcOn;
        int qmFormat ;
        Bool qmReplaceAll ;
        Bool qmMOPACAddConfigChrg;

        Bool qmBondOn;
        char qmBondColumn[16];
        Bool qmBondDist;
        int qmBondValType;
        char qmBondValueTypeS[16];
        char qmBondSchemeS[16] ;
        int qmBondScheme;
        Bool qmPCSwitchOn;
        char qmPCSwitchTypeS[16];
        int qmPCSwitchType;
        char qmPCSchemeS[16];
        int qmPCScheme;
        int qmSimsPerNode;

        Bool qmVDW ;
        Bool qmNoPC ;
        Bool qmElecEmbed ;
        int qmPCSelFreq ;
        Bool qmCustomPCSel;

        Bool qmLSSOn ;
        int qmLSSFreq ;
        char qmLSSResname[5] ;
        char qmLSSModeS[16];
        int qmLSSMode;

        Bool qmCSMD;
        char qmCSMDFile[NAMD_FILENAME_BUFFER_SIZE];

        int qmEnergyOutFreq ;
        int qmOutFreq ;
        int qmPosOutFreq ;

  Bool printBadContacts;        //print indices of bad contacts being moved downhill

  //gbis implicit solvent parameters
  Bool GBISOn;                    //do generalized born implicit solvent
  BigReal fsMax;
  Bool GBISserOn;                 //do generalized born implicit solvent serial
  BigReal solvent_dielectric;     //epsilon_s
  BigReal coulomb_radius_offset;  //rho_0
  BigReal kappa;      //debye screening length; k = sqrt(ion concentration mol/L ) / 0.304
  BigReal ion_concentration;
	BigReal gbis_delta;							//three parameters for born radius calc
	BigReal gbis_beta;
	BigReal gbis_gamma;
	BigReal alpha_cutoff;						//pairwise cutoff for integrating born radius
	BigReal alpha_max;								//maximum allowable born radius
  Bool LCPOOn;                    //do LCPO SASA for GBSA
  BigReal surface_tension;        //surface tension (kcal/mol/Ang^2) for LCPO

        Bool drudeOn;       // Perform integration of Drude oscillators?
        Bool drudeHardWallOn;  // Apply maximum Drude bond length restriction?
        BigReal drudeTemp;  // (low) temperature for freezing Drude oscillators
        BigReal drudeDamping;    // Langevin damping coefficient (1/ps)
                                 //   defaults to langevinDamping
        BigReal drudeBondLen;    // Length beyond which to apply quartic
                                 //   restraining potential to Drude bond
        BigReal drudeBondConst;  // Force constant for restraining potential
	BigReal drudeNbtholeCut;             // Radius of thole pair interaction

	Bool pairInteractionOn;		//  Calculate pair interactions?
	int pairInteractionGroup1;	//  Interaction group 1.
	int pairInteractionGroup2;	//  Interaction group 2.
        Bool pairInteractionSelf;       //  Compute just within group.

        Bool cosAngles;    // Can some angles be cos-based
	Bool globalForcesOn;		//  Are global forces present?
	Bool tclForcesOn;		//  Are Tcl forces present?
#ifdef NAMD_TCL
        Bool tclIsThreaded;             //  Is Tcl library thread-safe?
#endif
	Bool tclBCOn;			//  Are Tcl boundary forces present
	char *tclBCScript;		//  Script defining tclBC calcforces
	char tclBCArgs[NAMD_FILENAME_BUFFER_SIZE];		//  Extra args for calcforces command
	Bool freeEnergyOn;		//  Doing free energy perturbation?
	Bool miscForcesOn;		//  Using misc forces?
	Bool colvarsOn;         //  Using the colvars module?

	Bool fixedAtomsOn;		//  Are there fixed atoms?
	Bool fixedAtomsForces;		//  Calculate forces anyway?
	Bool fixedAtomsForceOutput; // Output fixed forces?

	Bool langevinOnAtStartup;	//  Ensure that langevin is set up properly
	Bool langevinOn;		//  Flag TRUE-> langevin dynamics active
	BigReal langevinTemp;		//  Temperature for Langevin dynamics
	BigReal langevinDamping;	//  Damping coefficient (1/ps)
	Bool langevinHydrogen;		//  Flag TRUE-> apply to hydrogens
	Bool langevin_useBAOAB;		//  Flag TRUE-> use the experimental BAOAB integrator for NVT instead of the BBK one
					//  See Leimkuhler and Matthews (AMRX 2012); implemented in NAMD by CM June2012

	// BEGIN LA
	Bool loweAndersenOn;		//  Flag TRUE-> Lowe-Andersen dynamics active
	BigReal loweAndersenTemp;	//  Temperature for Lowe-Andersen dynamics
	BigReal loweAndersenRate;	//  Collision frequency for Lowe-Andersen dynamics (1/ps)
	BigReal loweAndersenCutoff;	//  Cutoff radius for Lowe-Andersen dynamics
	// END LA

	Bool globalOn;			//  Flag TRUE-> use global integrator
	Bool dihedralOn;		//  Flag TRUE-> dihedral dynamics active
	Bool COLDOn;			//  Flag TRUE-> constrained overdamped
					//  langevin dynamics active
	BigReal COLDRate;		//  Damping coefficient for COLD.
	BigReal COLDTemp;		//  Temperature for COLD.

	Bool tCoupleOn;			//  Flag TRUE-> Temperature coupling
					//  active
	BigReal tCoupleTemp;		//  Temperature for temp coupling

  Bool stochRescaleOn;
  /**<  Flag TRUE enables the stochastic velocity rescaling method of
   * Bussi, Donadio, and Parrinello.  The method is an extension to the
   * Berendsen thermostat, where a stochastic update to the velocity
   * rescaling produces correct sampling of the NVT ensemble.
   */

  BigReal stochRescaleTemp;
  /**< Temperature for stochastic velocity rescaling. */

  BigReal stochRescalePeriod;
  /**< Timescale (ps) for stochastic velocity rescaling. */

  int stochRescaleFreq;
  /**< How frequently (time steps) stochastic velocity rescaling occurs. */

  Bool stochRescaleHeat;
  /**< Flag TRUE enables calculation and reporting of heat transfer and work.
    *  The computation is _cumulative_ from step = firstTimestep.
    */

	int rescaleFreq;		//  Velocity rescale frequency
	BigReal rescaleTemp;		//  Temperature to rescale to

        Bool accelMDOn;                 //  Perform accelerated MD
        Bool accelMDdihe;               //  Apply boost to the dihedral potential
        Bool accelMDdual;               //  dual boost mode
        Bool accelMDDebugOn;            //  Debugging accelerated MD
        BigReal accelMDFirstStep;       //  First aMD step
        BigReal accelMDLastStep;        //  Last aMD step
        int accelMDOutFreq;             //  aMD output frequency
        BigReal accelMDE;               //  aMD E
        BigReal accelMDalpha;           //  aMD alpha
        BigReal accelMDTE;              //  E for total potential in the dual boost mode
        BigReal accelMDTalpha;          //  alpha for total potential in the dual boost mode

	Bool accelMDG;                  //  Perform Gaussian accelMD calculation
	int accelMDGiE;                 //  Flag to set the mode iE in Gaussian accelMD
	int accelMDGcMDSteps;           //  Number of cMD steps
	int accelMDGEquiSteps;		//  Number of quilibration steps after adding boost potential
	int accelMDGcMDPrepSteps;	//  Number of preparation cMD steps
	int accelMDGEquiPrepSteps;	//  Number of preparation equilibration steps
        int accelMDGStatWindow;         //  Number of steps to calc avg and std
	BigReal accelMDGSigma0P;	//  upper limit of std of total potential
	BigReal accelMDGSigma0D;	//  upper limit of std of dihedral potential
	Bool accelMDGRestart;		//  Flag to set use restart file in Gaussian accelMD
  char accelMDGRestartFile[NAMD_FILENAME_BUFFER_SIZE];  //  restart file name
	Bool accelMDGresetVaftercmd;	//  Flag to reset potential after first accelMDGcMDSteps steps

        /* Begin Adaptive Temperature Sampling */
        Bool adaptTempOn;                      //  is adaptTempOn
        Bool adaptTempDebug;                   //  Debuggin adaptive temperature sampling
        int adaptTempFirstStep;                //  First adaptTemp step
        int adaptTempLastStep;                 //  Last adaptTemp step
        int adaptTempOutFreq;                  //  adaptTemp output frequency
        int adaptTempFreq;                     //  Steps between adaptTemp updates
        BigReal adaptTempTmin;                 //  Lower temperature bound
        BigReal adaptTempTmax;                 //  Upper temperature bound
        BigReal adaptTempAutoDt;               //  Auto jump size. Value determines upper bound, adaotTempDt determines lower bound
        int adaptTempBins;                     //  Number of bins to store average energy values
        BigReal adaptTempDt;                   //  timestep for adaptTemp updates - only affects Temperature random walk
        BigReal adaptTempCgamma;               //  Cgamma variable for adaptive bin averaging Cgamma = 0 is normal Averaging. 1 > Cgamma >= 0
        Bool adaptTempLangevin;                //  Couple to Langevin Thermostat
        Bool adaptTempRescale;                 //  Couple to Vel. Rescaling
        char adaptTempInFile[NAMD_FILENAME_BUFFER_SIZE];             //  Restart information for adaptTemp to read
        char adaptTempRestartFile[NAMD_FILENAME_BUFFER_SIZE];        //  File to write restart information
        int  adaptTempRestartFreq;             //  Frequency of writing restart output
        Bool adaptTempRandom;                  //  Do we assign random temperatures when we step out of [Tmin,Tmax]?
        /* End Adaptive Temperature Sampling */

	int reassignFreq;		//  Velocity reassignment frequency
	BigReal reassignTemp;		//  Temperature to reassign to
	BigReal reassignIncr;		//  Added to reassignTemp each time
	BigReal reassignHold;		//  Hold reassignTemp at this value

	Bool useGroupPressure;		//  Use group rather than atomic
					//  quantities for pressure calc

	Bool excludeFromPressure;	//  Flag TRUE-> some atoms not rescaled

	Bool useFlexibleCell;		//  Use anisotropic cell fluctuations
	Bool useConstantArea;		//  x,y dimensions fixed.
	Bool useConstantRatio;		//  x,y ratio fixed.

  Bool fixCellDims; // fix the cell dimensions
  Bool fixCellDimX;
  Bool fixCellDimY;
  Bool fixCellDimZ;

	Bool berendsenPressureOn;	//  Berendsen pressure bath
	BigReal berendsenPressureTarget;
	BigReal berendsenPressureCompressibility;
	BigReal berendsenPressureRelaxationTime;
	int berendsenPressureFreq;

	Bool langevinPistonOn;		//  Langevin piston pressure control
	Bool langevinPistonBarrier;	//  Turn off to extrapolate cell
	BigReal langevinPistonTarget;
	BigReal langevinPistonPeriod;
	BigReal langevinPistonDecay;
	BigReal langevinPistonTemp;

	Bool multigratorOn;     // Multigrator temperature and/or pressure control
	BigReal multigratorPressureTarget;
	BigReal multigratorPressureRelaxationTime;
	int multigratorPressureFreq;
	BigReal multigratorTemperatureTarget;
	BigReal multigratorTemperatureRelaxationTime;
	int multigratorTemperatureFreq;
	int multigratorNoseHooverChainLength;

        BigReal surfaceTensionTarget;

        Bool pressureProfileOn;         // Compute lateral pressure profile?
        int pressureProfileSlabs;       // Number of slabs
        int pressureProfileFreq;        // How often to store profile data
        int pressureProfileAtomTypes;
        Bool pressureProfileEwaldOn;    // Compute Ewald contribution?
        int pressureProfileEwaldX;
        int pressureProfileEwaldY;
        int pressureProfileEwaldZ;

	zVector strainRate;
	zVector strainRate2; // off diagonal elements (xy, xz, yz)

	unsigned int randomSeed;	//  Seed for random number generator

	Bool FMAOn;                     //  Flag TRUE-> FMA active
	int FMALevels;			//  Number of Levels for FMA
	int FMAMp;			//  Number of multipole terms for FMA
	Bool FMAFFTOn;			//  FFT on/off flag for FMA
	int FMAFFTBlock;		//  FFT blocking factor for FMA

	Bool fullDirectOn;		//  Should direct calculations of
					//  full electrostatics be performed?

        Bool MSMOn;      // enable MSM (multilevel summation method)
                         // for long-range electrostatics

        int MSMQuality;  // choose MSM quality 0 (low) - 3 (high), using
                         // optimal combination of approximation and splitting
                         // defaults to "low" for fastest performance

        int MSMApprox;   // choose MSM approximation
                         // defaults to "cubic" (low) for fastest performance

        int MSMSplit;    // choose MSM splitting function
                         // defaults to "Taylor2" (low) for fastest performance

        int MSMLevels;   // select number of MSM levels
                         // default (0) adapts number of levels to the
                         // system for fastest performance

        int MSMBlockSizeX;  // controls size of parallel work decomposition
        int MSMBlockSizeY;  // controls size of parallel work decomposition
        int MSMBlockSizeZ;  // controls size of parallel work decomposition

        BigReal MSMGridSpacing;  // defaults to 2.5 A, best for atomic systems

        BigReal MSMPadding;      // pad grid along non-periodic boundaries
                                 // defaults to 2.5 A
                                 // increase if atoms are drifting beyond
                                 // edge of grid, which will terminate
                                 // simulation prematurely

        BigReal MSMxmin;  // define extent of non-periodic boundaries
        BigReal MSMxmax;
        BigReal MSMymin;
        BigReal MSMymax;
        BigReal MSMzmin;
        BigReal MSMzmax;

        Bool MsmSerialOn;   // use serial MSM solver for testing

        Bool FMMOn;
        int FMMLevels;
        BigReal FMMPadding;

	Bool PMEOn;			//  Flag TRUE -> PME active
	BigReal PMETolerance;		//  Direct space tolerance
        BigReal PMEEwaldCoefficient;    //  From tolerance and cutoff
	int PMEInterpOrder;		//  Order of interpolation
	int PMEGridSizeX;		//  No. of grid points in x dim
	int PMEGridSizeY;		//  No. of grid points in y dim
	int PMEGridSizeZ;		//  No. of grid points in z dim
	BigReal PMEGridSpacing;		//  Maximum spacing between points
	int PMEProcessors;		//  No. of processors to use
	int PMEMinSlices;		//  Min slices per PME slab
	int PMEMinPoints;		//  Min points per PME pencil
	Bool PMEBarrier;		//  Use barrier before sendTrans
	int PMEPencils;			//  Size of pencil grid in each dim
	int PMEPencilsX;		//  Size of pencil grid in X dim
	int PMEPencilsY;		//  Size of pencil grid in Y dim
	int PMEPencilsZ;		//  Size of pencil grid in Z dim
	int PMEPencilsYLayout;		//  Y pencil layout strategy
	int PMEPencilsXLayout;		//  X pencil layout strategy
        int PMESendOrder;		//  Message ordering strategy
        Bool PMEOffload;		//  Offload reciprocal sum to accelerator

	Bool useDPME;			//  Flag TRUE -> old DPME code

	Bool usePMECUDA;
        /**< Flag TRUE to use the PME CUDA version.
         * Default is TRUE if running on 4 nodes or less. */

	Bool useCUDA2;
        /**< Flag TRUE to use the second generation nonbonded
         * CUDA kernels developed by Antti-Pekka Hynninen.
         * Default is TRUE. */

        int bondedCUDA;
        /**< Bitmask for calculating bonded interactions on GPU.
         * Default is 255, with the following bit position settings:
         * 1 -> bonds, 2 -> angles, 4 -> dihedrals, 8 -> impropers,
         * 16 -> exclusions, 32 -> crossterms. */

        Bool useCUDAdisable;
        /**< Flag TRUE to automatically disable individual CUDA kernels that
         * are incompatible with simulation options.  Default is TRUE.
         * Specifically, disable CUDA bonds for Drude oscillator simulation
         * and disable CUDA dihedral and CUDA crossterms kernels for
         * accelMDdihe and accelMDdual simulation.
         * Set FALSE to manually control kernel use for development. */

	Bool FFTWEstimate;
	Bool FFTWPatient;
	Bool FFTWUseWisdom;
	char FFTWWisdomFile[NAMD_FILENAME_BUFFER_SIZE];
	char *FFTWWisdomString;

        #ifdef OPENATOM_VERSION
        Bool openatom;                  // Flag TRUE -> OpenAtom QM/MM active
        #endif // OPENATOM_VERSION

	Bool minimizeCGOn;		//  Flag TRUE-> CG minimization active
        Bool minVerbose;                //  Flag TRUE-> print extra minimization data
	BigReal minTinyStep;		//  Minimization parameter
	BigReal minBabyStep;		//  Minimization parameter
	BigReal minLineGoal;		//  Minimization parameter
	Bool minimizeOn;		//  Flag TRUE-> minimization active
	BigReal maximumMove;		//  Maximum movement per timestep
					//  during minimization

	Bool sphericalBCOn;		//  Flag TRUE-> spherical boundary
					//  conditions are active
	zVector sphericalCenter;		//  Center specified by user
	BigReal sphericalBCk1;		//  First force constant for
					//  spherical BC
	BigReal sphericalBCk2;		//  Second force constant for
					//  spherical BC
	BigReal sphericalBCr1;		//  First radius for spherical BC
	BigReal sphericalBCr2;		//  Second radius for spherical BC
	int sphericalBCexp1;		//  First radius for spherical BC
	int sphericalBCexp2;		//  Second radius for spherical BC

        Bool cylindricalBCOn;           //  Flag TRUE->cylindrical boundary
                                        //  conditions are active
        zVector cylindricalCenter;
	char cylindricalBCAxis;		//  'x', 'y', or 'z'
        BigReal cylindricalBCr1;
        BigReal cylindricalBCr2;
        BigReal cylindricalBCl1;
        BigReal cylindricalBCl2;
        int cylindricalBCexp1;
        int cylindricalBCexp2;
        BigReal cylindricalBCk1;
        BigReal cylindricalBCk2;

	Bool eFieldOn;                  //  Should a electric field be applied
	Bool eFieldNormalized;          //  Is eField vector scaled by cell basis vectors
	zVector eField;                 //  Electric field vector to be applied
	BigReal eFieldFreq;		// Frequency of the electric field
	BigReal eFieldPhase;		// Phase phi, cos(w*t-phi*PI/180)

	Bool stirOn;                   // Should a stirring torque be applied
	char stirFilename[NAMD_FILENAME_BUFFER_SIZE];	       // Stirring filename (atoms marked)
	//do the below two even needed to be defined?
	BigReal stirStartingTheta;     // Stir starting theta offset
	BigReal stirVel;               // Stir angular velocity
	BigReal stirK;                 // Stir force harmonic spring constant
        zVector stirAxis;              // Direction of stir axis
        zVector stirPivot;             // Pivot point of stir axis

	Bool extraBondsOn;		// read extra bonded forces
	Bool extraBondsCosAngles;       // extra angles are cosine-based
	Bool extraBondsCosAnglesSetByUser; // did the user set this explicitly

	Bool consForceOn;		//  Should constant force be applied
  char consForceFile[NAMD_FILENAME_BUFFER_SIZE];
	BigReal consForceScaling;

	int outputEnergies;		//  Number of timesteps between energy
					//  outputs

	int outputMomenta;		//  Number of timesteps between momentum
					//  outputs

	int outputTiming;		//  Number of timesteps between timing
					//  outputs

	int outputCudaTiming;		//  Number of timesteps between timing
					//  outputs of CUDA code

	int outputPressure;		//  Number of timesteps between pressure
					//  tensor outputs

	Bool mergeCrossterms;		//  Merge crossterm energy w/ dihedrals

	int firstTimestep;		//  Starting timestep.  Will be 0 unless
					//  restarting a simulation

	MTSChoices MTSAlgorithm;	//  What multiple timestep algorithm
					//  to use

	int longSplitting;		//  What electrostatic splitting
					//  to use

	Bool ignoreMass;		//  Mass < 3.5 does not indicate hydrogen, etc.

	int splitPatch;			// How are patches determined?
	BigReal hgroupCutoff;		// what is the added hydrogen margin?

	int mollyOn;			// mollify long range forces?
	BigReal mollyTol;		// error tolerance for molly
	int mollyIter;			// max number of iterations for molly

        int rigidBonds;                 // what type of rigid bonds to hydrogens
                                        // none, all, or only water

        BigReal rigidTol;               // error tolerance for rigid bonds
        int rigidIter;                  // Number of NR iterations
	int rigidDie;			// die if rigidTol not achieved

	Bool useSettle;			// Use SETTLE; requires rigid waters

	Bool testOn;			//  Do tests rather than simulation
	Bool commOnly;			//  Don't do any force evaluations
	Bool statsOn;			//  Don't do any force evaluations

	int totalAtoms;			//  Total Number of atoms in simulation
        int maxSelfPart;                // maximum number of self partitions
                                        // that a patch can be split into
        int maxPairPart;                // maximum number of pair partitions
                                        // that a patch can be split into
        int numAtomsSelf;               // maximum number of atoms in a single
                                        // self-compute
        int numAtomsSelf2;              // maximum number of atoms in a pair compute
                                        // in the presence of twoAwayX,Y,Z options
        int numAtomsPair;               // maximum number of atoms in a single
                                        // pair-compute
        int numAtomsPair2;              // maximum number of atoms in a single
                                        // pair-compute
        int minAtomsPerPatch;           // minimum average atoms per patch
                                        //  (may create larger patches)
	int maxExclusionFlags;		// maximum size of exclusion check list
					// for any given atom
	Bool outputPatchDetails;	// print number of atoms per patch
        Bool staticAtomAssignment;      // never migrate atoms
        Bool replicaUniformPatchGrids;  // same patch grid size on all replicas

	//
        // hydrogen bond simulation parameters
        //

        // should the hydrogen bond term be used?  If FALSE, all other
	// hydrogen bond parameters are unnecessary in simulation.
	Bool HydrogenBonds;

	// should the antecedent atom be used in the calculation of hbonds?
	Bool useAntecedent;

	// exponents used in hydrogen bond energy function:
	//   aaAngleExp = exp for H-A-AA angle term (n)
	//   haAngleExp = exp for D-H-A angle term (m)
	//   distAttExp = exp for attractive A-D distance term (j)
	//   distRepExp = exp for repulsive A-D distance term (i)
	int aaAngleExp, haAngleExp, distAttExp, distRepExp;

	// cutoff D-H-A angle, and on/off angles for switch fcn (in degrees)
	BigReal dhaCutoffAngle, dhaOnAngle, dhaOffAngle;

	// cutoff distance for D-A separation in hbonds (in Angstroms), and
	// on/off distances for hbond radial term switching function
	BigReal daCutoffDist, daOnDist, daOffDist;

	// IMD parameters
	int IMDon;    // enable IMD
	int IMDport;  // port on which to listen for connections
 	int IMDfreq;  // frequency at which coordinates will be available
        int IMDwait;  // if true, pause the simulation when there is no
                      // connection
 	int IMDignore;  // IMD connection does not influence simulation
                        // only sends coordinates and energies to VMD
 	int IMDignoreForces;  // Only the Forces are ignored. Finish, Pause and Resume are enabled


        // AMBER options
        Bool amberOn; // FLAG TRUE-> amber force field is used
        Bool readExclusions; // FLAG TRUE-> Read exclusions from parm file
        BigReal vdwscale14; //  Scaling factor for 1-4 VDW interactions

	// GROMACS options
	Bool gromacsOn; // FLAG TRUE -> gromacs-style force field is used

	// OPLS options
	Bool vdwGeometricSigma;  // Lennard-J sigma uses geometric mean

	// ScriptTcl argument passing
	BigReal scriptArg1;
	BigReal scriptArg2;
	BigReal scriptArg3;
	BigReal scriptArg4;
	BigReal scriptArg5;
        int scriptIntArg1;
        int scriptIntArg2;
        char scriptStringArg1[128];
        char scriptStringArg2[128];

	Bool useCompressedPsf;
	Bool genCompressedPsf;

	Bool usePluginIO;

        Bool mallocTest;
        Bool printExclusions;

	//default value is -1
	int proxySendSpanningTree;
	int proxyRecvSpanningTree;

    int proxyTreeBranchFactor;


    //fields needed for Parallel IO Input
    int numinputprocs;
    char *binAtomFile;
    char *binCoorFile;
    char *binVelFile;
    char *binRefFile;

    //fields needed for Parallel IO Output
    int numoutputprocs;
    int numoutputwrts;

	char computeMapFilename[NAMD_FILENAME_BUFFER_SIZE];		//  store compute map
        Bool storeComputeMap;
        Bool loadComputeMap;

        // MIC-specific parameters
        int mic_hostSplit;
        int mic_numParts_self_p1;
        int mic_numParts_pair_p1;
        int mic_numParts_pair_p2;
        int mic_unloadMICPEs;
        int mic_deviceThreshold;
        int mic_singleKernel;

	// AVX-512 Tiles optimizations
	Bool useAVXTiles;

public:

        SimParameters() : mgridforcelist(), parseopts(0) {};
        SimParameters(ConfigList *c, char *&cwd) : mgridforcelist(), parseopts(0) {
	  initialize_config_data(c,cwd);
	};
        ~SimParameters() {};

	void initialize_config_data(ConfigList *, char *&cwd);
					//  Initialize SimParameters data
					//  from the ConfigList object
	void send_SimParameters(MOStream *);
					//  Used by the master process
					//  to send the paramters to
					//  the other processors
	void receive_SimParameters(MIStream *);
					//  Used by the other processors
					//  to receive the data from the
					//  master process
	void scriptSet(const char *, const char *);
					//  Set parameters at run time
	void close_dcdfile();  // *** implemented in Output.C ***
	void close_veldcdfile();  // *** implemented in Output.C ***
        static void nonbonded_select();
        static void pme_select();

	int isSendSpanningTreeOn(){ return proxySendSpanningTree == 1; }
	int isSendSpanningTreeUnset() { return proxySendSpanningTree == -1; }
	int isRecvSpanningTreeOn(){ return proxyRecvSpanningTree == 1; }
	int isRecvSpanningTreeUnset() { return proxyRecvSpanningTree == -1; }

        char* getfromparseopts(const char* name, char *outbuf);
        int istrueinparseopts(const char* name);
        int issetinparseopts(const char* name);

       	void readExtendedSystem(const char *filename, Lattice *latptr=0);
private:
        ParseOptions *parseopts;

	void config_parser(ParseOptions &opts);

	void config_parser_basic(ParseOptions &opts);
	void config_parser_fileio(ParseOptions &opts);
	void config_parser_fullelect(ParseOptions &opts);
	void config_parser_methods(ParseOptions &opts);
	void config_parser_constraints(ParseOptions &opts);
#ifdef OPENATOM_VERSION
	void config_parser_openatom(ParseOptions &opts);
#endif //OPENATOM_VERSION
	/* BEGIN gf */
	void config_parser_gridforce(ParseOptions &opts);
	/* END gf */
	void config_parser_movdrag(ParseOptions &opts);
	void config_parser_rotdrag(ParseOptions &opts);
	void config_parser_constorque(ParseOptions &opts);
	void config_parser_boundary(ParseOptions &opts);
	void config_parser_misc(ParseOptions &opts);
   	void config_parser_mgridforce(ParseOptions &opts);
     	void parse_mgrid_string_param(ConfigList *config,
     	                              const char *fieldname, char** dest);
     	void parse_mgrid_params(ConfigList *config);
       	void print_mgrid_params();

	void check_config(ParseOptions &opts, ConfigList *config, char *&cwd);

	void print_config(ParseOptions &opts, ConfigList *config, char *&cwd);

	void create_output_directories(const char *dirname);

	int fmaFrequency;		//  outdated parameter name
	char loadBalancer[64];		//  Load balancer
	char loadStrategy[64];		//  Load balancing strategy

};

#endif

