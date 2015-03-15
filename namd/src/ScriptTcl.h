/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Modifies SimParameters settings during run.
*/

#ifndef SCRIPTTCL_H
#define SCRIPTTCL_H

#include "converse.h"
#include "NamdTypes.h"
#include "Broadcasts.h"

#ifdef NAMD_TCL
#define USE_COMPAT_CONST
#include <tcl.h>
#endif

class ConfigList;
class NamdState;

class ScriptTcl {
public:
  ScriptTcl();
  ~ScriptTcl();
  void eval(char *script);
  void load(char *scriptFile);
#ifdef NAMD_TCL
  void run();
#else
  void run(char *scriptFile);
#endif
  void measure(Vector *);
private:
  char *scriptFile;
  ConfigList *config;
  NamdState *state;
  void suspend(void);
  int runWasCalled;
  int initWasCalled;
  void barrier();
  void initcheck();
  void reinitAtoms(const char *basename=0);
  SimpleBroadcastObject<int> scriptBarrier;
  int barrierStep;
  void runController(int task);
  void setParameter(const char* param, const char* value);
  void setParameter(const char* param, int value);
  friend class DataExchanger;
  int eval(const char *script, const char **resultPtr);
#ifdef NAMD_TCL
  friend class Controller;
  friend class GlobalMasterTcl;
  friend class colvarproxy_namd;
  Tcl_Interp *interp;
  static int Tcl_exit(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_abort(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_numPes(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_numNodes(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_numPhysicalNodes(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_numReplicas(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_myReplica(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_replicaEval(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_replicaYield(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_replicaSendrecv(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_replicaSend(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_replicaRecv(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_replicaBarrier(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_replicaAtomSendrecv(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_replicaAtomSend(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_replicaAtomRecv(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_stdout(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_print(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_config(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_isset_config(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_istrue_config(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_param(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_isset_param(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_istrue_param(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_reinitvels(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_rescalevels(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_run(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_minimize(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_move(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_moveallby(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_output(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_measure(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_colvarbias(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_colvarvalue(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_colvars(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_colvarfreq(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_checkpoint(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_revert(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_checkpointReplica(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_callback(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_reinitatoms(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_coorfile(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_dumpbench(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_consForceConfig(ClientData, Tcl_Interp *, int, Tcl_Obj *const objv[]);
  static int Tcl_reloadCharges(ClientData, Tcl_Interp *, int, char **);
  static int Tcl_reloadGridforceGrid(ClientData, Tcl_Interp *, int, char **);	// BEGIN, END gf
  char *callbackname;
  void doCallback(const char *labels, const char *data);
  int doCallback() { return ! ! callbackname; }
  char *measure_command;
  int measure_result;
#endif
};

#endif

