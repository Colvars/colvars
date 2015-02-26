/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Modifies SimParameters settings during run.
*/

#include "InfoStream.h"
#include "BackEnd.h"
#include "ScriptTcl.h"
#include "Broadcasts.h"
#include "ConfigList.h"
#include "Node.h"
#include "PDB.h"
#include "WorkDistrib.h"
#include "NamdState.h"
#include "Controller.h"
#include "SimParameters.h"
#include "Thread.h"
#include "ProcessorPrivate.h"
#include "PatchMgr.h"
#include "PatchMap.h"
#include "Measure.h"
#include "colvarmodule.h"
#include "colvarscript.h"
#include "DumpBench.h"
#include <errno.h>
#include <stdio.h>
#include <ctype.h>  // for isspace
#ifndef WIN32
#include <strings.h>
#endif

#include "qd.h"

#ifdef NAMD_TCL
#define USE_COMPAT_CONST
#include <tcl.h>
#endif
#include "TclCommands.h"

#include "ProcessorPrivate.h"
#include "DataExchanger.h"

//#define DEBUGM
#define MIN_DEBUG_LEVEL 4
#include "Debug.h"

#include "molfile_plugin.h"
#include "libmolfile_plugin.h"

static molfile_plugin_t *dcdplugin;
static int register_cb(void *v, vmdplugin_t *p) {
	dcdplugin = (molfile_plugin_t *)p;
	return 0;
}

//
// XXX static and global variables are unsafe for shared memory builds.
//
static int numatoms;
static void *filehandle;
static float *coords;
static Vector *vcoords;


void ScriptTcl::suspend() {
  BackEnd::suspend();
}

void ScriptTcl::barrier() {
  BackEnd::barrier();
}

void ScriptTcl::initcheck() {
  if ( initWasCalled == 0 ) {
#ifdef NAMD_TCL
    CkPrintf("TCL: Suspending until startup complete.\n");
    Tcl_CreateCommand(interp, "param", Tcl_param,
      (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "unknown", Tcl_param,
      (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "isset", Tcl_isset_param,
      (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
    Tcl_CreateCommand(interp, "istrue", Tcl_istrue_param,
      (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
#endif
    initWasCalled = 1;

    state->configListInit(config);
    Node::Object()->saveMolDataPointers(state);
#ifdef NAMD_TCL
    SimParameters *simParams = Node::Object()->simParameters;
    simParams->tclIsThreaded =
      ! ! Tcl_GetVar2(interp, "tcl_platform", "threaded", TCL_GLOBAL_ONLY);
#endif
    Node::messageStartUp();
    suspend();
  }
}

void ScriptTcl::runController(int task) {
  scriptBarrier.publish(barrierStep++,task);
  suspend();
#ifdef NAMD_TCL
  if ( task == SCRIPT_RUN || task == SCRIPT_CONTINUE || task == SCRIPT_MINIMIZE  ) {
    doCallback(state->callback_labelstring.c_str(),
               state->callback_valuestring.c_str());
  }
#endif
}

void ScriptTcl::setParameter(const char* param, const char* value) {
  ScriptParamMsg *msg = new ScriptParamMsg;
  strncpy(msg->param,param,MAX_SCRIPT_PARAM_SIZE);
  strncpy(msg->value,value,MAX_SCRIPT_PARAM_SIZE);
  (CProxy_Node(CkpvAccess(BOCclass_group).node)).scriptParam(msg);
  barrier();
}

void ScriptTcl::setParameter(const char* param, int value) {
  ScriptParamMsg *msg = new ScriptParamMsg;
  strncpy(msg->param,param,MAX_SCRIPT_PARAM_SIZE);
  sprintf(msg->value,"%d",value);
  (CProxy_Node(CkpvAccess(BOCclass_group).node)).scriptParam(msg);
  barrier();
}

void ScriptTcl::reinitAtoms(const char *basename) {
  Node::Object()->workDistrib->reinitAtoms(basename);
  barrier();
}

#ifdef NAMD_TCL

int ScriptTcl::Tcl_exit(ClientData clientData,
	Tcl_Interp *, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  CkPrintf("TCL: Exiting due to exit command.\n");
  script->runController(SCRIPT_END);
  BackEnd::exit();
  return TCL_OK;
}

int ScriptTcl::Tcl_abort(ClientData,
	Tcl_Interp *, int argc, char *argv[]) {
  Tcl_DString msg;
  Tcl_DStringInit(&msg);
  Tcl_DStringAppend(&msg,"TCL:",-1);
  for ( int i = 1; i < argc; ++i ) {
    Tcl_DStringAppend(&msg," ",-1);
    Tcl_DStringAppend(&msg,argv[i],-1);
  }
  NAMD_die(Tcl_DStringValue(&msg));
  Tcl_DStringFree(&msg);
  return TCL_OK;
}

int ScriptTcl::Tcl_numPes(ClientData, Tcl_Interp *interp, int argc, char **) {
  if ( argc > 1 ) {
    Tcl_SetResult(interp,"no arguments needed",TCL_VOLATILE);
    return TCL_ERROR;
  }
  Tcl_SetObjResult(interp, Tcl_NewIntObj(CkNumPes()));
  return TCL_OK;
}

int ScriptTcl::Tcl_numNodes(ClientData, Tcl_Interp *interp, int argc, char **) {
  if ( argc > 1 ) {
    Tcl_SetResult(interp,"no arguments needed",TCL_VOLATILE);
    return TCL_ERROR;
  }
  Tcl_SetObjResult(interp, Tcl_NewIntObj(CkNumNodes()));
  return TCL_OK;
}

int ScriptTcl::Tcl_numPhysicalNodes(ClientData, Tcl_Interp *interp, int argc, char **) {
  if ( argc > 1 ) {
    Tcl_SetResult(interp,"no arguments needed",TCL_VOLATILE);
    return TCL_ERROR;
  }
  Tcl_SetObjResult(interp, Tcl_NewIntObj(CmiNumPhysicalNodes()));
  return TCL_OK;
}

int ScriptTcl::Tcl_numReplicas(ClientData, Tcl_Interp *interp, int argc, char **) {
  if ( argc > 1 ) {
    Tcl_SetResult(interp,"no arguments needed",TCL_VOLATILE);
    return TCL_ERROR;
  }
  Tcl_SetObjResult(interp, Tcl_NewIntObj(CmiNumPartitions()));
  return TCL_OK;
}

int ScriptTcl::Tcl_myReplica(ClientData, Tcl_Interp *interp, int argc, char **) {
  if ( argc > 1 ) {
    Tcl_SetResult(interp,"no arguments needed",TCL_VOLATILE);
    return TCL_ERROR;
  }
  Tcl_SetObjResult(interp, Tcl_NewIntObj(CmiMyPartition()));
  return TCL_OK;
}

#define CHECK_REPLICA(REP) do {\
  if ( (REP) < 0 ) { \
    Tcl_SetResult(interp,"negative replica index",TCL_VOLATILE); \
    return TCL_ERROR; \
  } \
  if ( (REP) >= CmiNumPartitions() ) { \
    Tcl_SetResult(interp,"non-existent replica index",TCL_VOLATILE); \
    return TCL_ERROR; \
  } \
} while ( 0 )

int ScriptTcl::Tcl_replicaEval(ClientData, Tcl_Interp *interp, int argc, char **argv) {
  if ( argc != 3 ) {
    Tcl_SetResult(interp,"args: dest script",TCL_VOLATILE);
    return TCL_ERROR;
  }
  int dest = atoi(argv[1]);
  CHECK_REPLICA(dest);
#if CMK_HAS_PARTITION
  Tcl_DString recvstr;
  Tcl_DStringInit(&recvstr);
  DataMessage *recvMsg = NULL;
  replica_eval(argv[2], dest, CkMyPe(), &recvMsg);
  CmiAssert(recvMsg != NULL);
  int code = recvMsg->code;
  Tcl_DStringAppend(&recvstr, recvMsg->data, recvMsg->size);
  Tcl_DStringResult(interp, &recvstr);
  Tcl_DStringFree(&recvstr);
  CmiFree(recvMsg);
  return code;
#else
  return Tcl_EvalEx(interp,argv[2],-1,TCL_EVAL_GLOBAL);
#endif
}


int ScriptTcl::Tcl_replicaSendrecv(ClientData, Tcl_Interp *interp, int argc, char **argv) {
  if ( argc < 3 || argc > 4 ) {
    Tcl_SetResult(interp,"args: data dest ?source?",TCL_VOLATILE);
    return TCL_ERROR;
  }
  Tcl_DString recvstr;
  Tcl_DStringInit(&recvstr);
  int sendcount = strlen(argv[1]);
  int recvcount = 0;
  int dest = atoi(argv[2]);
  int source = -1;
  if ( argc > 3 ) source = atoi(argv[3]);
#if CMK_HAS_PARTITION
  if(dest == CmiMyPartition()) {
    Tcl_DStringSetLength(&recvstr,sendcount);
    memcpy(Tcl_DStringValue(&recvstr),argv[1],sendcount);
  } else {
    DataMessage *recvMsg = NULL;
    replica_sendRecv(argv[1], sendcount, dest, CkMyPe(), &recvMsg, source, CkMyPe());
    CmiAssert(recvMsg != NULL);
    Tcl_DStringAppend(&recvstr, recvMsg->data, recvMsg->size);
    CmiFree(recvMsg);
  }
#endif
  Tcl_DStringResult(interp, &recvstr);
  Tcl_DStringFree(&recvstr);
  return TCL_OK;
}

int ScriptTcl::Tcl_replicaSend(ClientData, Tcl_Interp *interp, int argc, char **argv) {
  if ( argc != 3 ) {
    Tcl_SetResult(interp,"args: data dest",TCL_VOLATILE);
    return TCL_ERROR;
  }
  int sendcount = strlen(argv[1]);
  int dest = atoi(argv[2]);
#if CMK_HAS_PARTITION
  replica_send(argv[1], sendcount, dest, CkMyPe());
#endif
  return TCL_OK;
}

int ScriptTcl::Tcl_replicaRecv(ClientData, Tcl_Interp *interp, int argc, char **argv) {
  if (argc != 2 ) {
    Tcl_SetResult(interp,"args: source",TCL_VOLATILE);
    return TCL_ERROR;
  }
  Tcl_DString recvstr;
  Tcl_DStringInit(&recvstr);
  int recvcount = 0;
  int source = atoi(argv[1]);
#if CMK_HAS_PARTITION
  DataMessage *recvMsg = NULL;
  replica_recv(&recvMsg, source, CkMyPe());
  CmiAssert(recvMsg != NULL);
  Tcl_DStringAppend(&recvstr, recvMsg->data, recvMsg->size);
  CmiFree(recvMsg);
#endif
  Tcl_DStringResult(interp, &recvstr);
  Tcl_DStringFree(&recvstr);
  return TCL_OK;
}

int ScriptTcl::Tcl_replicaBarrier(ClientData, Tcl_Interp *interp, int argc, char **) {
  if ( argc > 1 ) {
    Tcl_SetResult(interp,"no arguments needed",TCL_VOLATILE);
    return TCL_ERROR;
  }
#if CMK_HAS_PARTITION
  replica_barrier();
#endif
  return TCL_OK;
}

int ScriptTcl::Tcl_replicaAtomSendrecv(ClientData clientData, Tcl_Interp *interp, int argc, char **argv) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  if ( ! Node::Object()->simParameters->replicaUniformPatchGrids ) {
    Tcl_SetResult(interp,"replicaUniformPatchGrids is required for atom exchange",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( argc < 2 || argc > 3 ) {
    Tcl_SetResult(interp,"bad arg count; args: dest ?source?",TCL_VOLATILE);
    return TCL_ERROR;
  }
  int dest = -1;
  if ( sscanf(argv[1], "%d", &dest) != 1 ) {
    Tcl_SetResult(interp,"bad dest; args: dest ?source?",TCL_VOLATILE);
    return TCL_ERROR;
  }
  int source = -1;
  if ( argc == 3 ) {
    if ( sscanf(argv[2], "%d", &source) != 1 ) {
      Tcl_SetResult(interp,"bad source; args: dest ?source?",TCL_VOLATILE);
      return TCL_ERROR;
    }
  }

#if CMK_HAS_PARTITION
  if(dest != CmiMyPartition()) {
    DataMessage *recvMsg = NULL;
    replica_sendRecv((char*)&(script->state->lattice), sizeof(Lattice), dest, CkMyPe(), &recvMsg, source, CkMyPe());
    CmiAssert(recvMsg != NULL);
    memcpy(&(script->state->lattice), recvMsg->data, recvMsg->size);
    CmiFree(recvMsg);
  }
#endif

  char str[40];
  sprintf(str, "%d", dest);
  script->setParameter("scriptArg1", str);
  sprintf(str, "%d", source);
  script->setParameter("scriptArg2", str);

  CkpvAccess(_qd)->create(2 * PatchMap::Object()->numPatches());

  script->runController(SCRIPT_ATOMSENDRECV);

#if CMK_HAS_PARTITION
  if(dest != CmiMyPartition()) {
    DataMessage *recvMsg = NULL;
    ControllerState *cstate = script->state->controller;
    replica_sendRecv((char*)cstate, sizeof(ControllerState), dest, CkMyPe(), &recvMsg, source, CkMyPe());
    CmiAssert(recvMsg != NULL);
    memcpy(cstate, recvMsg->data, recvMsg->size);
    CmiFree(recvMsg);
  }
#endif

  return TCL_OK;
}

int ScriptTcl::Tcl_replicaAtomSend(ClientData clientData, Tcl_Interp *interp, int argc, char **argv) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  if ( ! Node::Object()->simParameters->replicaUniformPatchGrids ) {
    Tcl_SetResult(interp,"replicaUniformPatchGrids is required for atom exchange",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( argc != 2 ) {
    Tcl_SetResult(interp,"bad arg count; args: dest",TCL_VOLATILE);
    return TCL_ERROR;
  }
  int dest = -1;
  if ( sscanf(argv[1], "%d", &dest) != 1 ) {
    Tcl_SetResult(interp,"bad dest; args: dest",TCL_VOLATILE);
    return TCL_ERROR;
  }

#if CMK_HAS_PARTITION
  replica_send((char*)&(script->state->lattice), sizeof(Lattice), dest, CkMyPe());
#endif

  char str[40];
  sprintf(str, "%d", dest);
  script->setParameter("scriptArg1", str);

  CkpvAccess(_qd)->create(PatchMap::Object()->numPatches());

  script->runController(SCRIPT_ATOMSEND);

#if CMK_HAS_PARTITION
  ControllerState *cstate = script->state->controller;
  replica_send((char*)cstate, sizeof(ControllerState), dest, CkMyPe());
#endif

  return TCL_OK;
}

int ScriptTcl::Tcl_replicaAtomRecv(ClientData clientData, Tcl_Interp *interp, int argc, char **argv) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  if ( ! Node::Object()->simParameters->replicaUniformPatchGrids ) {
    Tcl_SetResult(interp,"replicaUniformPatchGrids is required for atom exchange",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( argc > 2 ) {
    Tcl_SetResult(interp,"bad arg count; args: ?source?",TCL_VOLATILE);
    return TCL_ERROR;
  }
  int source = -1;
  if ( argc == 2 ) {
    if ( sscanf(argv[1], "%d", &source) != 1 ) {
      Tcl_SetResult(interp,"bad source; args: ?source?",TCL_VOLATILE);
      return TCL_ERROR;
    }
  }

#if CMK_HAS_PARTITION
  DataMessage *recvMsg = NULL;
  replica_recv(&recvMsg, source, CkMyPe());
  CmiAssert(recvMsg != NULL);
  memcpy(&(script->state->lattice), recvMsg->data, recvMsg->size);
  CmiFree(recvMsg);
#endif

  char str[40];
  sprintf(str, "%d", source);
  script->setParameter("scriptArg2", str);

  CkpvAccess(_qd)->create(PatchMap::Object()->numPatches());

  script->runController(SCRIPT_ATOMRECV);

#if CMK_HAS_PARTITION
  recvMsg = NULL;
  ControllerState *cstate = script->state->controller;
  replica_recv(&recvMsg, source, CkMyPe());
  CmiAssert(recvMsg != NULL);
  memcpy(cstate, recvMsg->data, recvMsg->size);
  CmiFree(recvMsg);
#endif

  return TCL_OK;
}


int ScriptTcl::Tcl_stdout(ClientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  if (argc != 2) {
    Tcl_SetResult(interp, "wrong # args", TCL_VOLATILE);
    return TCL_ERROR;
  }

  char *filename= argv[1];
  CkPrintf("TCL: redirecting stdout to file %s\n", filename);

  if ( ! freopen(filename, "a", stdout) ) {
    Tcl_SetResult(interp, strerror(errno), TCL_VOLATILE);
    return TCL_ERROR;
  }
  return TCL_OK;
}

int ScriptTcl::Tcl_print(ClientData,
	Tcl_Interp *, int argc, char *argv[]) {
  Tcl_DString msg;
  Tcl_DStringInit(&msg);
  for ( int i = 1; i < argc; ++i ) {
    Tcl_DStringAppend(&msg," ",-1);
    Tcl_DStringAppend(&msg,argv[i],-1);
  }
  CkPrintf("TCL:%s\n",Tcl_DStringValue(&msg));
  Tcl_DStringFree(&msg);
  return TCL_OK;
}

int ScriptTcl::Tcl_config(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {

// Needs to handle the following cases as passed in by Tcl:
//    name data #comment
//    name=data #comment
//    name= data #comment
//    name =data #comment
//    name = data #comment
//    name data1 data2 data3 #comment
//    name=data1 data2 data3 #comment
//    name= data1 data2 data3 #comment
//    name =data1 data2 data3 #comment
//    name = data1 data2 data3 #comment
//    name { data1 data2 data3 } #comment
//    name { data1 data2 data3 } #comment
//    name { data1 data2 # data3 } #comment
//    name {data1 data2 # data3 } #comment
// Do not try to handle "data#comment" in any form.
// The '#' start of any comments will *always* be a new argv.
// The name will *always* be contained in argv[1].

  // allocate storage for data string
  int arglen = 1;  int ai;
  for (ai=1; ai<argc; ++ai) { arglen += strlen(argv[ai]) + 1; }
  char *data = new char[arglen];  *data = 0;

  // find the end of the name
  char *name, *s;
  name = argv[1];
  for ( s = name; *s && *s != '='; ++s );

  // eliminate any comment
  for (ai=2; ai<argc; ++ai) { if (argv[ai][0] == '#') argc = ai; }

  // concatenate all the data items
  ai = 2;
  if ( *s ) { *s = 0; ++s; strcat(data,s); }  // name=data or name=
  else if ( ai < argc && argv[ai][0] == '=' ) {  // name =data or name =
    strcat(data,argv[ai]+1);
    ++ai;
  }
  for ( ; ai<argc; ++ai) {
    if ( data[0] ) { strcat(data," "); }
    strcat(data,argv[ai]);
  }

  if ( ! *name ) {
    delete [] data;
    Tcl_SetResult(interp,"error parsing config file",TCL_VOLATILE);
    return TCL_ERROR;
  }

  ScriptTcl *script = (ScriptTcl *)clientData;

  if ( *data ) {
    script->config->add_element( name, strlen(name), data, strlen(data) );
    delete [] data;
    return TCL_OK;
  }
  delete [] data;

  StringList *strlist = script->config->find(name);
  if ( ! strlist ) {
    Tcl_SetResult(interp,"error parsing config file",TCL_VOLATILE);
    return TCL_ERROR;
  }
  Tcl_SetResult(interp,strlist->data,TCL_VOLATILE);
  return TCL_OK;
}

int ScriptTcl::Tcl_isset_config(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  if (argc != 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }

  char *param = argv[1];
  ScriptTcl *script = (ScriptTcl *)clientData;
  StringList *strlist = script->config->find(param);
  Tcl_SetResult(interp, (char*)(strlist ? "1" : "0"), TCL_VOLATILE);
  return TCL_OK;
}

static int atoBool(const char *s)
{
   if (!strcasecmp(s, "on")) return 1;
   if (!strcasecmp(s, "off")) return 0;
   if (!strcasecmp(s, "true")) return 1;
   if (!strcasecmp(s, "false")) return 0;
   if (!strcasecmp(s, "yes")) return 1;
   if (!strcasecmp(s, "no")) return 0;
   if (!strcasecmp(s, "1")) return 1;
   if (!strcasecmp(s, "0")) return 0;
   return -1;
}

int ScriptTcl::Tcl_istrue_config(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  if (argc != 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }

  char *param = argv[1];
  ScriptTcl *script = (ScriptTcl *)clientData;
  StringList *strlist = script->config->find(param);
  if ( ! strlist ) {
    Tcl_SetResult(interp,"parameter value is not set",TCL_VOLATILE);
    return TCL_ERROR;
  }
  int val = atoBool(strlist->data);
  if ( val < 0 ) {
    Tcl_SetResult(interp,"parameter value is not boolean",TCL_VOLATILE);
    return TCL_ERROR;
  }
  Tcl_SetResult(interp, (char*)(val ? "1" : "0"), TCL_VOLATILE);
  return TCL_OK;
}

int ScriptTcl::Tcl_istrue_param(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  if (argc != 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }

  char *param = argv[1];
  SimParameters *simParams = Node::Object()->simParameters;
  int val = simParams->istrueinparseopts(param);
  if ( val == -1 ) {
    Tcl_SetResult(interp,"unknown parameter",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( val == -2 ) {
    Tcl_SetResult(interp,"parameter is not boolean",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( val == -3 ) {
    Tcl_SetResult(interp,"parameter value is not set",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( val != 0 && val != 1 ) {
    Tcl_SetResult(interp,"bug in Tcl_istrue_param",TCL_VOLATILE);
    return TCL_ERROR;
  }
  Tcl_SetResult(interp, (char*)(val ? "1" : "0"), TCL_VOLATILE);
  return TCL_OK;
}

int ScriptTcl::Tcl_isset_param(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  if (argc != 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }

  char *param = argv[1];
  SimParameters *simParams = Node::Object()->simParameters;
  int val = simParams->issetinparseopts(param);
  if ( val < 0 ) {
    Tcl_SetResult(interp,"unknown parameter",TCL_VOLATILE);
    return TCL_ERROR;
  }
  Tcl_SetResult(interp, (char*)(val ? "1" : "0"), TCL_VOLATILE);
  return TCL_OK;
}

int ScriptTcl::Tcl_param(ClientData clientData,
        Tcl_Interp *interp, int argc, char *argv[]) {
  if (argc != 2 && argc != 3 && argc != 5) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }

  char *param = argv[1];
  if ( strlen(param) + 1 > MAX_SCRIPT_PARAM_SIZE ) {
    Tcl_SetResult(interp,"parameter name too long",TCL_VOLATILE);
    return TCL_ERROR;
  }

  if ( argc == 2 ) { // get param value
    char buf[MAX_SCRIPT_PARAM_SIZE];
    SimParameters *simParams = Node::Object()->simParameters;
    char *result = simParams->getfromparseopts(param,buf);
    if ( result ) {
      Tcl_SetResult(interp, result,TCL_VOLATILE);
      return TCL_OK;
    } else {
      Tcl_SetResult(interp,"unknown parameter",TCL_VOLATILE);
      return TCL_ERROR;
    }
  }

  char value[MAX_SCRIPT_PARAM_SIZE];
  int arglen = strlen(argv[2]) + 1;
  if ( argc == 5 ) arglen += strlen(argv[3]) + strlen(argv[4]) + 2;
  if ( arglen > MAX_SCRIPT_PARAM_SIZE ) {
    Tcl_SetResult(interp,"parameter value too long",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( argc == 3 ) sprintf(value,"%s",argv[2]);
  if ( argc == 5 ) sprintf(value,"%s %s %s",argv[2],argv[3],argv[4]);

  iout << "TCL: Setting parameter " << param << " to " << value << "\n" << endi;

  ScriptTcl *script = (ScriptTcl *)clientData;
  script->setParameter(param,value);

  return TCL_OK;
}

int ScriptTcl::Tcl_reinitvels(ClientData clientData,
        Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  if (argc != 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  char *temp = argv[1];

  script->setParameter("initialTemp",temp);

  script->runController(SCRIPT_REINITVELS);

  return TCL_OK;
}

int ScriptTcl::Tcl_rescalevels(ClientData clientData,
        Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  if (argc != 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  char *factor = argv[1];

  script->setParameter("scriptArg1",factor);

  script->runController(SCRIPT_RESCALEVELS);

  return TCL_OK;
}

int ScriptTcl::Tcl_run(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  if (argc < 2) {
    Tcl_SetResult(interp,"too few args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if (argc > 3) {
    Tcl_SetResult(interp,"too many args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  int norepeat = 0;
  if (argc == 3) {
    if ( ! strcmp(argv[1], "norepeat") ) {
      if ( script->runWasCalled ) { norepeat = 1; }
    } else {
      Tcl_SetResult(interp,
        "first arg not norepeat",TCL_VOLATILE);
      return TCL_ERROR;
    }
  }
  int numstepsarg = argc-1;
  int numsteps;
  if (Tcl_GetInt(interp,argv[numstepsarg],&numsteps) != TCL_OK) {
    return TCL_ERROR;
  }
  if (numsteps < 0) {
    Tcl_SetResult(interp,"number of steps must be non-negative",TCL_VOLATILE);
    return TCL_ERROR;
  }
  SimParameters *simParams = Node::Object()->simParameters;
  if (numsteps && simParams->firstTimestep % simParams->stepsPerCycle) {
    Tcl_SetResult(interp,"firstTimestep must be a multiple of stepsPerCycle",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if (numsteps % simParams->stepsPerCycle) {
    Tcl_SetResult(interp,"number of steps must be a multiple of stepsPerCycle",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( simParams->N != simParams->firstTimestep ) {
    iout << "TCL: Original numsteps " << simParams->N
         << " will be ignored.\n";
  }
  iout << "TCL: Running for " << numsteps << " steps";
  if ( norepeat ) iout << " without repeating first step";
  iout << "\n" << endi;

  script->setParameter("numsteps",simParams->firstTimestep + numsteps);

  script->runController(norepeat ? SCRIPT_CONTINUE : SCRIPT_RUN);
  script->runWasCalled = 1;

  script->setParameter("firsttimestep",simParams->N);

  return TCL_OK;
}

int ScriptTcl::Tcl_minimize(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  if (argc != 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  int numsteps;
  if (Tcl_GetInt(interp,argv[1],&numsteps) != TCL_OK) {
    return TCL_ERROR;
  }
  if (numsteps < 0) {
    Tcl_SetResult(interp,"number of steps must be non-negative",TCL_VOLATILE);
    return TCL_ERROR;
  }
  SimParameters *simParams = Node::Object()->simParameters;
  if (numsteps && simParams->firstTimestep % simParams->stepsPerCycle) {
    Tcl_SetResult(interp,"firstTimestep must be a multiple of stepsPerCycle",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if (numsteps % simParams->stepsPerCycle) {
    Tcl_SetResult(interp,"number of steps must be a multiple of stepsPerCycle",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( simParams->N != simParams->firstTimestep ) {
    iout << "TCL: Original numsteps " << simParams->N
         << " will be ignored.\n";
  }
  iout << "TCL: Minimizing for " << numsteps << " steps\n" << endi;

  script->setParameter("numsteps",simParams->firstTimestep + numsteps);

  script->runController(SCRIPT_MINIMIZE);
  script->runWasCalled = 1;

  script->setParameter("firsttimestep",simParams->N);

  return TCL_OK;
}

// move all atoms by a given vector
int ScriptTcl::Tcl_moveallby(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  if (argc != 2) {
    Tcl_SetResult(interp, "wrong # args", TCL_VOLATILE);
    return TCL_ERROR;
  }
  char **fstring;
  int fnum;
  double x, y, z;
  if (Tcl_SplitList(interp, argv[1], &fnum, &fstring) != TCL_OK)
    return TCL_ERROR;
  if ( (fnum != 3) ||
       (Tcl_GetDouble(interp, fstring[0],&x) != TCL_OK) ||
       (Tcl_GetDouble(interp, fstring[1],&y) != TCL_OK) ||
       (Tcl_GetDouble(interp, fstring[2],&z) != TCL_OK) ) {
    Tcl_SetResult(interp,"argument not a vector",TCL_VOLATILE);
    Tcl_Free((char*)fstring);
    return TCL_ERROR;
  }
  Tcl_Free((char*)fstring);

  MoveAllByMsg *msg = new MoveAllByMsg;
  msg->offset = Vector(x,y,z);
  (CProxy_PatchMgr(CkpvAccess(BOCclass_group).patchMgr)).moveAllBy(msg);

  script->barrier();
  return TCL_OK;
}

int ScriptTcl::Tcl_move(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  if (argc != 4) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  char **fstring;  int fnum;  int atomid;  int moveto;  double x, y, z;
  if (Tcl_GetInt(interp,argv[1],&atomid) != TCL_OK) return TCL_ERROR;
  if (argv[2][0]=='t' && argv[2][1]=='o' && argv[2][2]==0) moveto = 1;
  else if (argv[2][0]=='b' && argv[2][1]=='y' && argv[2][2]==0) moveto = 0;
  else {
    Tcl_SetResult(interp,"syntax is 'move <id> to|by {<x> <y> <z>}'",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if (Tcl_SplitList(interp, argv[3], &fnum, &fstring) != TCL_OK) {
    return TCL_ERROR;
  }
  if ( (fnum != 3) ||
       (Tcl_GetDouble(interp, fstring[0],&x) != TCL_OK) ||
       (Tcl_GetDouble(interp, fstring[1],&y) != TCL_OK) ||
       (Tcl_GetDouble(interp, fstring[2],&z) != TCL_OK) ) {
    Tcl_SetResult(interp,"third argument not a vector",TCL_VOLATILE);
    Tcl_Free((char*)fstring);
    return TCL_ERROR;
  }
  Tcl_Free((char*)fstring);

  SimParameters *simParams = Node::Object()->simParameters;

  iout << "TCL: Moving atom " << atomid << " ";
  if ( moveto ) iout << "to"; else iout << "by";
  iout << " " << Vector(x,y,z) << ".\n" << endi;

  MoveAtomMsg *msg = new MoveAtomMsg;
  msg->atomid = atomid - 1;
  msg->moveto = moveto;
  msg->coord = Vector(x,y,z);
  (CProxy_PatchMgr(CkpvAccess(BOCclass_group).patchMgr)).moveAtom(msg);

  script->barrier();

  return TCL_OK;
}

int ScriptTcl::Tcl_output(ClientData clientData,
        Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  if (argc < 2) {
    Tcl_SetResult(interp,"too few args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if (argc > 3) {
    Tcl_SetResult(interp,"too many args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  int filenamearg = argc-1;
  if (strlen(argv[filenamearg]) > MAX_SCRIPT_PARAM_SIZE) {
    Tcl_SetResult(interp,"file name too long",TCL_VOLATILE);
    return TCL_ERROR;
  }
  int dorestart = 1;
  int doforces = 0;
  if (argc == 3) {
    if ( ! strcmp(argv[1], "withforces") ) {
      doforces = 1;
    } else if ( ! strcmp(argv[1], "onlyforces") ) {
      dorestart = 0;
      doforces = 1;
    } else {
      Tcl_SetResult(interp,
        "first arg not withforces or onlyforces",TCL_VOLATILE);
      return TCL_ERROR;
    }
  }

  SimParameters *simParams = Node::Object()->simParameters;

  char oldname[MAX_SCRIPT_PARAM_SIZE+1];
  strncpy(oldname,simParams->outputFilename,MAX_SCRIPT_PARAM_SIZE);

  script->setParameter("outputname",argv[filenamearg]);

  iout << "TCL: Writing to files with basename " <<
		simParams->outputFilename << ".\n" << endi;

  if ( doforces && ! script->runWasCalled ) NAMD_die(
    "No forces to output; must call run or minimize first.");

  if ( dorestart ) script->runController(SCRIPT_OUTPUT);
  if ( doforces ) script->runController(SCRIPT_FORCEOUTPUT);

  script->setParameter("outputname",oldname);

  return TCL_OK;
}

void ScriptTcl::measure(Vector *c) {
  Measure::createCommands(interp);
  Node::Object()->coords = c;
  measure_result = Tcl_Eval(interp,measure_command);
  Node::Object()->coords = 0;
  Measure::deleteCommands(interp);
}

int ScriptTcl::Tcl_measure(ClientData clientData,
        Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  if (argc != 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  script->measure_command = argv[1];

  script->runController(SCRIPT_MEASURE);

  return script->measure_result;
}

// NOTE: This interface is DEPRECATED
// Please use the "cv bias" interface instead:

// Replace "colvarbias changeconfig" with:
// cv bias <name> delete
// cv config <new_config_string>

// Replace "colvarbias energydiff" with:
// cv bias config <config_string_with_tempBias>
// set ediff [expr [cv bias tempBias energy] - [cv bias refBias energy]]
// cv bias tempBias delete

int ScriptTcl::Tcl_colvarbias(ClientData clientData,
        Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  if (argc < 4 || argc % 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  colvarmodule *colvars = Node::Object()->colvars;
  if ( ! colvars ) {
    Tcl_SetResult(interp,"colvars module not active",TCL_VOLATILE);
    return TCL_ERROR;
  }
  if ( ! strcmp(argv[1],"changeconfig") ) {
    for ( int i=2; i<argc; i+=2 ) {
      std::string name(argv[i]);
      std::string conf(argv[i+1]);
      colvars->change_configuration(name,conf);
    }
    return TCL_OK;
  } else if ( ! strcmp(argv[1],"energydiff") ) {
    if ( ! script->runWasCalled ) {
      Tcl_SetResult(interp,"energydiff requires a previous timestep",TCL_VOLATILE);
      return TCL_ERROR;
    }
    double ediff = 0.;
    for ( int i=2; i<argc; i+=2 ) {
      std::string name(argv[i]);
      std::string conf(argv[i+1]);
      ediff += colvars->energy_difference(name,conf);
    }
    Tcl_SetObjResult(interp, Tcl_NewDoubleObj(ediff));
    return TCL_OK;
  } else {
    Tcl_SetResult(interp,"unknown colvarbias operation",TCL_VOLATILE);
    return TCL_ERROR;
  }
}

// NOTE: This interface is DEPRECATED
// Please use the "cv colvar" interface instead

int ScriptTcl::Tcl_colvarvalue(ClientData clientData,
        Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  if (argc != 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  colvarmodule *colvars = Node::Object()->colvars;
  if ( ! colvars ) {
    Tcl_SetResult(interp,"colvars module not active",TCL_VOLATILE);
    return TCL_ERROR;
  }
  // Pass the colvarvalue to Tcl
  std::string name(argv[1]);
  std::string value = colvars->read_colvar(name);
  // Process from a colvar list to a Tcl compatible list
  size_t found;
  do {
    found = value.find("(");
    if (found != std::string::npos) {
      value.replace(found, 1, " ");
    } else {
      break;
    }
  } while (true);
  do {
    found = value.find(")");
    if (found != std::string::npos) {
      value.replace(found, 1, " ");
    } else {
      break;
    }
  } while (true);
  do {
    found = value.find(",");
    if (found != std::string::npos) {
      value.replace(found, 1, " ");
    } else {
      break;
    }
  } while (true);
  // Send the result to Tcl
  Tcl_DString recvstr;
  Tcl_DStringInit(&recvstr);
  Tcl_DStringAppend(&recvstr,value.c_str(), -1);
  Tcl_DStringResult(interp, &recvstr);
  Tcl_DStringFree(&recvstr);
  return TCL_OK;
}

int ScriptTcl::Tcl_colvarfreq(ClientData clientData,
        Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  if (argc != 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }
  colvarmodule *colvars = Node::Object()->colvars;
  if ( ! colvars ) {
    Tcl_SetResult(interp,"colvars module not active",TCL_VOLATILE);
    return TCL_ERROR;
  }
  int new_freq;
  if (Tcl_GetInt(interp,argv[1],&new_freq) != TCL_OK) {
    return TCL_ERROR;
  }
  colvars->cv_traj_freq = new_freq;
  return TCL_OK;
}

int ScriptTcl::Tcl_colvars (ClientData clientData,
        Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  colvarmodule *colvars = Node::Object()->colvars;
  if ( ! colvars ) {
    Tcl_SetResult(interp,"colvars module not active",TCL_VOLATILE);
    return TCL_ERROR;
  }
  int retval = colvars->proxy->script->run(argc, (char const **) argv);
  // use Tcl dynamic allocation to prevent having to copy the buffer
  // *twice* just because Tcl is missing const qualifiers for strings
  char *buf = Tcl_Alloc(colvars->proxy->script->result.length() + 1);
  strncpy (buf, colvars->proxy->script->result.c_str(), colvars->proxy->script->result.length() + 1);
  Tcl_SetResult(interp, buf, TCL_DYNAMIC);
  // Note: sometimes Tcl 8.5 will segfault here
  // (only on error conditions, apparently)
  // http://sourceforge.net/p/tcl/bugs/4677/
  // Fixed in Tcl 8.6

  if (retval == COLVARSCRIPT_OK && !cvm::get_error())
    return TCL_OK;
  else
    return TCL_ERROR;
}

int ScriptTcl::Tcl_checkpoint(ClientData clientData,
        Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  if (argc != 1) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }

  script->runController(SCRIPT_CHECKPOINT);

  return TCL_OK;
}

int ScriptTcl::Tcl_revert(ClientData clientData,
        Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  if (argc != 1) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }

  script->runController(SCRIPT_REVERT);

  return TCL_OK;
}

int ScriptTcl::Tcl_callback(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  if (argc != 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }

  delete [] script->callbackname;
  script->callbackname = new char[strlen(argv[1])+1];
  strcpy(script->callbackname,argv[1]);

  iout << "TCL: Reduction callback proc set to " <<
			script->callbackname << "\n" << endi;

  return TCL_OK;
}

void ScriptTcl::doCallback(const char *labels, const char *data) {
  if ( ! callbackname ) return;
  int len = strlen(callbackname) + strlen(labels) + strlen(data) + 7;
  char *cmd = new char[len];
  sprintf(cmd, "%s {%s} {%s}", callbackname, labels, data);
  int rval = Tcl_Eval(interp,cmd);
  delete [] cmd;
  if (rval != TCL_OK) {
    const char *errorInfo = Tcl_GetVar(interp,"errorInfo",0);
    NAMD_die(errorInfo);
  }
}

extern void read_binary_coors(char *fname, PDB *pdbobj);

int ScriptTcl::Tcl_reinitatoms(ClientData clientData,
        Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  if (argc > 2) {
    Tcl_SetResult(interp,"wrong # args",TCL_VOLATILE);
    return TCL_ERROR;
  }

  if (argc == 1 ) {
    iout << "TCL: Reinitializing atom data\n" << endi;
    SimParameters *simParams = Node::Object()->simParameters;
    Controller *c = script->state->controller;
    script->state->lattice = c->origLattice;
    c->langevinPiston_strainRate = c->langevinPiston_origStrainRate;
    c->rescaleVelocities_sumTemps = 0;  c->rescaleVelocities_numTemps = 0;
    c->berendsenPressure_avg = 0; c->berendsenPressure_count = 0;
    SetLatticeMsg *msg = new SetLatticeMsg;
    msg->lattice = script->state->lattice;
    (CProxy_PatchMgr(CkpvAccess(BOCclass_group).patchMgr)).setLattice(msg);
    script->barrier();
    if ( ! simParams->binaryOutput ) {  // output may have overwritten data in PDB
      StringList *coordinateFilename = script->state->configList->find("bincoordinates");
      if ( coordinateFilename ) {
        read_binary_coors(coordinateFilename->data, script->state->pdb);
      } else if (coordinateFilename = script->state->configList->find("coordinates")) {
        PDB coordpdb(coordinateFilename->data);
        if ( coordpdb.num_atoms() != script->state->pdb->num_atoms() ) {
          NAMD_die("inconsistent atom count on re-reading coordinates pdb file");
        }
        Vector *positions = new Position[coordpdb.num_atoms()];
        coordpdb.get_all_positions(positions);
        script->state->pdb->set_all_positions(positions);
        delete [] positions;
      } else {
        iout << iWARN << "reinitatoms may fail if pdb-format output has occurred\n" << endi;
      }
    }
    script->reinitAtoms();
    return TCL_OK;
  }

  iout << "TCL: Reinitializing atom data from files with basename " << argv[1] << "\n" << endi;
  SimParameters *simParams = Node::Object()->simParameters;
  simParams->readExtendedSystem((std::string(argv[1])+".xsc").c_str(), &(script->state->lattice));
  Controller *c = script->state->controller;
  c->langevinPiston_strainRate =
      Tensor::symmetric(simParams->strainRate,simParams->strainRate2);
  c->rescaleVelocities_sumTemps = 0;  c->rescaleVelocities_numTemps = 0;
  c->berendsenPressure_avg = 0; c->berendsenPressure_count = 0;
  SetLatticeMsg *msg = new SetLatticeMsg;
  msg->lattice = script->state->lattice;
  (CProxy_PatchMgr(CkpvAccess(BOCclass_group).patchMgr)).setLattice(msg);
  script->barrier();
  script->reinitAtoms(argv[1]);

  return TCL_OK;
}

#define DEG2RAD 3.14159625359/180.0
#define UNITCELLSLOP 0.0001

static int get_lattice_from_ts(Lattice *lattice, const molfile_timestep_t *ts)
{
  // Check if valid unit cell data is contained in the timestep.  We don't
  // have any formalized way of doing this yet; for now, just check that
  // the length of the vector is greater than 1.
  if (ts->A <= 1 || ts->B <= 1 || ts->C <= 1) return 0;

  // convert from degrees to radians
  // Try to get exact results when the angles are exactly 90.
  double epsalpha = DEG2RAD*(ts->alpha-90.0);
  double epsbeta  = DEG2RAD*(ts->beta-90.0);
  double epsgamma = DEG2RAD*(ts->gamma-90.0);
  double cosAB = -sin(epsgamma);
  double sinAB = cos(epsgamma);
  double cosAC = -sin(epsbeta);
  double cosBC = -sin(epsalpha);

  // A will lie along the positive x axis.
  // B will lie in the x-y plane
  // The origin will be (0,0,0).
  Vector A(0), B(0), vecC(0);
  A.x = ts->A;
  B.x = ts->B*cosAB;
  B.y = ts->B*sinAB;
  //if (fabs(B.x) < UNITCELLSLOP) B.x = 0;
  //if (fabs(B.y) < UNITCELLSLOP) B.y = 0;
  vecC.x = ts->C * cosAC;
  vecC.y = (ts->B*ts->C*cosBC - B.x*vecC.x)/B.y;
  vecC.z = sqrt(ts->C*ts->C - vecC.x*vecC.x - vecC.y*vecC.y);
  //if (fabs(vecC.x) < UNITCELLSLOP) vecC.x = 0;
  //if (fabs(vecC.y) < UNITCELLSLOP) vecC.y = 0;
  //if (fabs(vecC.z) < UNITCELLSLOP) vecC.z = 0;
  lattice->set(A, B, vecC, Vector(0));
  return 1;
}

int ScriptTcl::Tcl_coorfile(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  if (argc == 4 && !strcmp(argv[1], "open")) {
    if (strcmp(argv[2], "dcd")) {
      NAMD_die("Sorry, coorfile presently supports only DCD files");
    }
    filehandle = dcdplugin->open_file_read(argv[3], "dcd", &numatoms);
    if (!filehandle) {
      Tcl_AppendResult(interp, "coorfile: Error opening file ", argv[3], NULL);
      return TCL_ERROR;
    }
    if (numatoms != Node::Object()->pdb->num_atoms()) {
      Tcl_AppendResult(interp, "Coordinate file ", argv[3],
        "\ncontains the wrong number of atoms.", NULL);
      return TCL_ERROR;
    }
    coords = new float[3*numatoms];
    vcoords = new Vector[3*numatoms];
    iout << iINFO << "Coordinate file " << argv[3] << " opened for reading.\n"
         << endi;
  } else if (argc == 2 && !strcmp(argv[1], "read")) {
    if (filehandle == NULL) {
      Tcl_AppendResult(interp, "coorfile read: Error, no file open for reading",
	NULL);
      return TCL_ERROR;
    }
    molfile_timestep_t ts;
    ts.coords = coords;
    int rc = dcdplugin->read_next_timestep(filehandle, numatoms, &ts);
    if (rc) {  // EOF
      Tcl_SetObjResult(interp, Tcl_NewIntObj(-1));
      return TCL_OK;
    }
    iout << iINFO << "Reading timestep from file.\n" << endi;
    Lattice lattice;
    if (get_lattice_from_ts(&lattice, &ts)) {
      iout << iINFO << "Updating unit cell from timestep.\n" << endi;
      if ( lattice.a_p() && ! script->state->lattice.a_p() ||
           lattice.b_p() && ! script->state->lattice.b_p() ||
           lattice.c_p() && ! script->state->lattice.c_p() ) {
        iout << iWARN << "Cell basis vectors should be specified before reading trajectory.\n" << endi;
      }
      // update Controller's lattice, but don't change the origin!
      Vector a(0.);  if ( script->state->lattice.a_p() ) a = lattice.a();
      Vector b(0.);  if ( script->state->lattice.b_p() ) b = lattice.b();
      Vector c(0.);  if ( script->state->lattice.c_p() ) c = lattice.c();
      script->state->lattice.set(a,b,c);
      SetLatticeMsg *msg = new SetLatticeMsg;
      msg->lattice = script->state->lattice;
      (CProxy_PatchMgr(CkpvAccess(BOCclass_group).patchMgr)).setLattice(msg);
      script->barrier();
    }
    for (int i=0; i<numatoms; i++) {
      vcoords[i].x = coords[3*i+0];
      vcoords[i].y = coords[3*i+1];
      vcoords[i].z = coords[3*i+2];
    }
    Node::Object()->pdb->set_all_positions(vcoords);
    script->reinitAtoms();
    Tcl_SetObjResult(interp, Tcl_NewIntObj(0));
  } else if (argc == 2 && !strcmp(argv[1], "close")) {
    if (!filehandle) {
      Tcl_AppendResult(interp, "coorfile close: No file opened for reading!",
        NULL);
      return TCL_OK;
    }
    iout << iINFO << "Closing coordinate file.\n" << endi;
    dcdplugin->close_file_read(filehandle);
    filehandle = NULL;
    delete [] coords;
    delete [] vcoords;

  } else if (argc ==2 && !strcmp(argv[1], "skip")) {
    if (filehandle == NULL) {
      Tcl_AppendResult(interp, "coorfile skip: Error, no file open for reading",
	NULL);
      return TCL_ERROR;
    }
    int rc = dcdplugin->read_next_timestep(filehandle, numatoms, NULL);
    if (rc) {  // EOF
      Tcl_SetObjResult(interp, Tcl_NewIntObj(-1));
      return TCL_OK;
    }
    Tcl_SetObjResult(interp, Tcl_NewIntObj(0));

  } else {
    NAMD_die("Unknown option passed to coorfile");
  }
  return TCL_OK;
}

int ScriptTcl::Tcl_dumpbench(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  if (argc != 2) {
    Tcl_AppendResult(interp, "usage: dumpbench <filename>", NULL);
    return TCL_ERROR;
  }

  if ( CkNumPes() != 1 ) {
    Tcl_AppendResult(interp, "multiple processors detected; dumpbench only works on serial runs", NULL);
    return TCL_ERROR;
  }

  FILE *file = fopen(argv[1],"w");
  if ( ! file ) {
    Tcl_AppendResult(interp, "dumpbench: error opening file ", argv[1], NULL);
    return TCL_ERROR;
  }

  if ( dumpbench(file) ) {
    Tcl_AppendResult(interp, "dumpbench: error dumping benchmark data", NULL);
    return TCL_ERROR;
  }

  fclose(file);

  Tcl_AppendResult(interp, "benchmark data written to file ", argv[1], NULL);
  return TCL_OK;
}

#include "ComputeConsForceMsgs.h"
// consforceconfig <atomids> <forces>
int ScriptTcl::Tcl_consForceConfig(ClientData clientData,
    Tcl_Interp *interp, int objc, Tcl_Obj *const objv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  if ( ! Node::Object()->simParameters->consForceOn ) {
    Tcl_AppendResult(interp, "consForceConfig requires constantForce on", NULL);
    return TCL_ERROR;
  }
  if (objc != 3) {
    Tcl_WrongNumArgs(interp, 1, objv, (char *)"<atomids> <forces>");
    return TCL_ERROR;
  }
  int natoms, nforces;
  Tcl_Obj **atomobjlist, **forceobjlist;
  if (Tcl_ListObjGetElements(interp, objv[1], &natoms, &atomobjlist) != TCL_OK ||
      Tcl_ListObjGetElements(interp, objv[2], &nforces, &forceobjlist) != TCL_OK) {
    return TCL_ERROR;
  }
  if (natoms != nforces) {
    Tcl_AppendResult(interp, (char *)"consforceconfig: atom list and force list not the same size!", NULL);
    return TCL_ERROR;
  }
  ComputeConsForceMsg *msg = new ComputeConsForceMsg;
  for (int i=0; i<natoms; i++) {
    int atomid;
    int nelem;
    Tcl_Obj **elemlist;
    Vector force;
    if (Tcl_GetIntFromObj(interp, atomobjlist[i], &atomid) != TCL_OK)
      return TCL_ERROR;
    if (Tcl_ListObjGetElements(interp, forceobjlist[i], &nelem, &elemlist) != TCL_OK)
      return TCL_ERROR;
    if (nelem != 3) {
      Tcl_AppendResult(interp, (char *)"consforceconfig: forces must have three elements", NULL);
      return TCL_ERROR;
    }
    if (Tcl_GetDoubleFromObj(interp, elemlist[0], &force.x) != TCL_OK ||
        Tcl_GetDoubleFromObj(interp, elemlist[1], &force.y) != TCL_OK ||
        Tcl_GetDoubleFromObj(interp, elemlist[2], &force.z) != TCL_OK) {
      return TCL_ERROR;
    }
    msg->aid.add(atomid);
    msg->f.add(force);
  }
  (CProxy_ComputeMgr(CkpvAccess(BOCclass_group).computeMgr)).recvComputeConsForceMsg(msg);
  return TCL_OK;
}

int ScriptTcl::Tcl_reloadCharges(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();
  if (argc != 2) {
    Tcl_AppendResult(interp, "usage: reloadCharges <filename>", NULL);
    return TCL_ERROR;
  }

  Node::Object()->reloadCharges(argv[1]);

  script->runController(SCRIPT_RELOADCHARGES);

  return TCL_OK;
}

// BEGIN gf
int ScriptTcl::Tcl_reloadGridforceGrid(ClientData clientData,
	Tcl_Interp *interp, int argc, char *argv[]) {
  ScriptTcl *script = (ScriptTcl *)clientData;
  script->initcheck();

  char *key = NULL;
  if (argc == 1) {
      // nothing ... key is NULL, then Node::reloadGridforceGrid uses the
      // default key, which is used internally when the gridforce*
      // keywords are used (as opposed to the mgridforce* keywords)
  } else if (argc == 2) {
      key = argv[1];
  } else {
      Tcl_AppendResult(interp, "usage: reloadGridforceGrid [<gridkey>]", NULL);
      return TCL_ERROR;
  }

  //(CProxy_Node(CkpvAccess(BOCclass_group).node)).reloadGridforceGrid(key);
  Node::Object()->reloadGridforceGrid(key);
  script->barrier();

  return TCL_OK;
}
// END gf

#endif  // NAMD_TCL


ScriptTcl::ScriptTcl() : scriptBarrier(scriptBarrierTag) {
  DebugM(3,"Constructing ScriptTcl\n");
#ifdef NAMD_TCL
  interp = 0;
  callbackname = 0;
#endif
  state = new NamdState;
  barrierStep = 0;

  molfile_dcdplugin_init();
  molfile_dcdplugin_register(NULL, register_cb);

  initWasCalled = 0;
  runWasCalled = 0;

#ifdef NAMD_TCL
  config = new ConfigList;

  // Create interpreter
  interp = Tcl_CreateInterp();
  tcl_vector_math_init(interp);
  Tcl_CreateCommand(interp, "exit", Tcl_exit,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "abort", Tcl_abort,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "numPes", Tcl_numPes,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "numNodes", Tcl_numNodes,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "numPhysicalNodes", Tcl_numPhysicalNodes,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "numReplicas", Tcl_numReplicas,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "myReplica", Tcl_myReplica,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "replicaEval", Tcl_replicaEval,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "replicaSendrecv", Tcl_replicaSendrecv,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "replicaSend", Tcl_replicaSend,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "replicaRecv", Tcl_replicaRecv,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "replicaBarrier", Tcl_replicaBarrier,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "replicaAtomSendrecv", Tcl_replicaAtomSendrecv,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "replicaAtomSend", Tcl_replicaAtomSend,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "replicaAtomRecv", Tcl_replicaAtomRecv,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "stdout", Tcl_stdout,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "print", Tcl_print,
    (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "unknown", Tcl_config,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "param", Tcl_config,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "isset", Tcl_isset_config,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "istrue", Tcl_istrue_config,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "run", Tcl_run,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "minimize", Tcl_minimize,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "move", Tcl_move,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "moveallby", Tcl_moveallby,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "output", Tcl_output,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "measure", Tcl_measure,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "colvarbias", Tcl_colvarbias,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "colvarvalue", Tcl_colvarvalue,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "cv", Tcl_colvars,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "colvarfreq", Tcl_colvarfreq,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "checkpoint", Tcl_checkpoint,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "revert", Tcl_revert,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "reinitvels", Tcl_reinitvels,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "rescalevels", Tcl_rescalevels,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "reinitatoms", Tcl_reinitatoms,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "callback", Tcl_callback,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "coorfile", Tcl_coorfile,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "dumpbench", Tcl_dumpbench,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateObjCommand(interp, "consForceConfig", Tcl_consForceConfig,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  Tcl_CreateCommand(interp, "reloadCharges", Tcl_reloadCharges,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  // BEGIN gf
  Tcl_CreateCommand(interp, "reloadGridforceGrid", Tcl_reloadGridforceGrid,
    (ClientData) this, (Tcl_CmdDeleteProc *) NULL);
  // END gf
#endif

}

int ScriptTcl::eval(const char *script, const char **resultPtr) {

#ifdef NAMD_TCL
  int code = Tcl_EvalEx(interp,script,-1,TCL_EVAL_GLOBAL);
  *resultPtr = Tcl_GetStringResult(interp);
  return code;
#else
  NAMD_bug("ScriptTcl::eval called without Tcl.");
#endif
}

void ScriptTcl::eval(char *script) {

#ifdef NAMD_TCL
  int code = Tcl_Eval(interp,script);
  const char *result = Tcl_GetStringResult(interp);
  if (*result != 0) CkPrintf("TCL: %s\n",result);
  if (code != TCL_OK) {
    const char *errorInfo = Tcl_GetVar(interp,"errorInfo",0);
    NAMD_die(errorInfo);
  }
#else
  NAMD_bug("ScriptTcl::eval called without Tcl.");
#endif

}

void ScriptTcl::load(char *scriptFile) {

#ifdef NAMD_TCL
  int code = Tcl_EvalFile(interp,scriptFile);
  const char *result = Tcl_GetStringResult(interp);
  if (*result != 0) CkPrintf("TCL: %s\n",result);
  if (code != TCL_OK) {
    const char *errorInfo = Tcl_GetVar(interp,"errorInfo",0);
    NAMD_die(errorInfo);
  }
#else
  NAMD_bug("ScriptTcl::load called without Tcl.");
#endif

}

#ifdef NAMD_TCL
void ScriptTcl::run() {
#else
void ScriptTcl::run(char *scriptFile) {

  if ( NULL == scriptFile || NULL == (config = new ConfigList(scriptFile)) ) {
    NAMD_die("Simulation config file is empty.");
  }
#endif

  if (runWasCalled == 0) {
    initcheck();
    SimParameters *simParams = Node::Object()->simParameters;
    if ( simParams->minimizeCGOn ) runController(SCRIPT_MINIMIZE);
    else runController(SCRIPT_RUN);
    runWasCalled = 1;
  }

  runController(SCRIPT_END);

}

ScriptTcl::~ScriptTcl() {
  DebugM(3,"Destructing ScriptTcl\n");
#ifdef NAMD_TCL
  if ( interp ) Tcl_DeleteInterp(interp);
  delete [] callbackname;
#endif

  molfile_dcdplugin_fini();
}

