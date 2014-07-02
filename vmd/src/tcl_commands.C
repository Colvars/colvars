/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2011 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: tcl_commands.C,v $
 *      $Author: johns $        $Locker:  $             $State: Exp $
 *      $Revision: 1.42 $       $Date: 2014/03/18 21:28:53 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *   Fundamental VMD Tcl text commands
 ***************************************************************************/

#include <stdlib.h>
#include <string.h>
#include "tcl.h"
#include "TclCommands.h"
#include "tcl_commands.h"
#include "config.h"
#include "utilities.h"
#include "CUDAKernels.h"
#include "WKFThreads.h"
#include "vmd.h"
#if defined(VMDCOLVARS)
#include "colvarproxy_vmd.h"
#endif

class VMDApp;

#define SIMPLE_TCL_OPT(string,result)       \
if (!strcmp(argv[1], string)) {             \
  Tcl_AppendResult(interp, result, NULL);   \
  return TCL_OK;                            \
}

static int vmdinfo_tcl(ClientData, Tcl_Interp *interp,
                       int argc, const char *argv[]) {
  VMDApp *app = (VMDApp *)Tcl_GetAssocData(interp, "VMDApp", NULL);

  if (argc == 2) {
    SIMPLE_TCL_OPT("version", VMDVERSION);
    SIMPLE_TCL_OPT("versionmsg", VERSION_MSG);
    SIMPLE_TCL_OPT("authors", VMD_AUTHORS);
    SIMPLE_TCL_OPT("arch", VMD_ARCH);
    SIMPLE_TCL_OPT("options", VMD_OPTIONS);
    SIMPLE_TCL_OPT("www", VMD_HOMEPAGE);
    SIMPLE_TCL_OPT("wwwhelp", VMD_HELPPAGE);


    // return the estimated amount of available physical memory
    if (!strcmp(argv[1], "freemem")) {
      long vmdcorefree = vmd_get_avail_physmem_mb();
      Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
      Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewIntObj(vmdcorefree));
      Tcl_SetObjResult(interp, tcl_result);
      return TCL_OK;
    }


    // return the number of available CPU cores
    if (!strcmp(argv[1], "numcpus")) {
#if defined(VMDTHREADS)
      int numcpus = wkf_thread_numprocessors();
#else
      int numcpus = 1;
#endif
      Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
      Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewIntObj(numcpus));
      Tcl_SetObjResult(interp, tcl_result);
      return TCL_OK;
    }


    // return the CPU affinity list for the VMD process
    if (!strcmp(argv[1], "cpuaffinity")) {
      int numcpus = -1;
      int *cpuaffinitylist = NULL;

#if defined(VMDTHREADS)
      cpuaffinitylist = wkf_cpu_affinitylist(&numcpus);
#endif
      if (numcpus > 0 && cpuaffinitylist != NULL) {
        int i;
        Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
        for (i=0; i<numcpus; i++)
          Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewIntObj(cpuaffinitylist[i]));
        Tcl_SetObjResult(interp, tcl_result);
        return TCL_OK;
      }

      if (cpuaffinitylist != NULL)
        free(cpuaffinitylist);

      Tcl_AppendResult(interp, "CPU affinity query unavailable on this platform", NULL);
      return TCL_ERROR;
    }


    // return the number of available CUDA devices
    if (!strcmp(argv[1], "numcudadevices")) {
      int numdevices;
#if defined(VMDCUDA)
      vmd_cuda_num_devices(&numdevices);
#else
      numdevices = 0;
#endif
      Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
      Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewIntObj(numdevices));
      Tcl_SetObjResult(interp, tcl_result);
      return TCL_OK;
    }

    // return the active display device (e.g. "text", "win", "cave", ...)
    if (!strcmp(argv[1], "dispdev")) {
      const char *disp = VMDgetDisplayTypeName();
      Tcl_AppendResult(interp, disp, NULL);
      return TCL_OK;
    }

    // return the MPI node name 
    if (!strcmp(argv[1], "nodename")) {
      Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
      Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewStringObj(app->par_name(), strlen(app->par_name())));
      Tcl_SetObjResult(interp, tcl_result);
      return TCL_OK;
    }  

    // return the MPI node rank 
    if (!strcmp(argv[1], "noderank")) {
      Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
      Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewIntObj(app->par_rank()));
      Tcl_SetObjResult(interp, tcl_result);
      return TCL_OK;
    }  

    // return the MPI node count
    if (!strcmp(argv[1], "nodecount")) {
      Tcl_Obj *tcl_result = Tcl_NewListObj(0, NULL);
      Tcl_ListObjAppendElement(interp, tcl_result, Tcl_NewIntObj(app->par_size()));
      Tcl_SetObjResult(interp, tcl_result);
      return TCL_OK;
    }  
  }

  Tcl_AppendResult(interp,
    "vmdinfo: version | versionmsg | authors | arch | \n"
    "freemem | numcpus | cpuaffinity | numcudadevices | \n"
    "displaytype | nodename | noderank | nodecount | \n"
    "options | www | wwwhelp", NULL);
  return TCL_ERROR;
}


int Vmd_Init(Tcl_Interp *interp) {
  VMDApp *app = (VMDApp *)Tcl_GetAssocData(interp, "VMDApp", NULL);

  Tcl_CreateCommand(interp,  "vmdinfo", vmdinfo_tcl,
        (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateCommand(interp,  "vmdbench", text_cmd_vmdbench,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateCommand(interp, "animate", text_cmd_animate,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateCommand(interp, "color", text_cmd_color,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateCommand(interp, "axes", text_cmd_axes,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateCommand(interp, "display", text_cmd_display,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateCommand(interp, "imd", text_cmd_imd,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateCommand(interp, "vmdcollab", text_cmd_collab,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateCommand(interp, "vmd_label", text_cmd_label,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateCommand(interp, "light", text_cmd_light,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateCommand(interp, "pointlight", text_cmd_point_light,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateCommand(interp, "material", text_cmd_material,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateCommand(interp, "vmd_menu", text_cmd_menu,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);
  
  Tcl_CreateCommand(interp, "stage", text_cmd_stage,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateCommand(interp, "light", text_cmd_light,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateCommand(interp, "user", text_cmd_user,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateCommand(interp, "mol", text_cmd_mol,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateCommand(interp, "molecule", text_cmd_mol,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateCommand(interp, "mouse", text_cmd_mouse,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateCommand(interp, "mobile", text_cmd_mobile,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateCommand(interp, "spaceball", text_cmd_spaceball,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateCommand(interp, "plugin", text_cmd_plugin,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateCommand(interp, "render", text_cmd_render,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

#if defined(VMDTK) && !defined(_MSC_VER)
  Tcl_CreateCommand(interp, "tkrender", text_cmd_tkrender,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);
#endif

  Tcl_CreateCommand(interp, "rock", text_cmd_rock,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateCommand(interp, "rotate", text_cmd_rotate,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateCommand(interp, "rotmat", text_cmd_rotmat,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateCommand(interp, "vmd_scale", text_cmd_scale,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateCommand(interp, "translate", text_cmd_translate,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateCommand(interp, "sleep", text_cmd_sleep,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

#if 1
  Tcl_CreateObjCommand(interp, "mdffi", obj_mdff_cc,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);
#endif
  
#if 0
  Tcl_CreateObjCommand(interp, "volgradient", obj_volgradient,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);
#endif

  Tcl_CreateCommand(interp, "tool", text_cmd_tool,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateObjCommand(interp,  "measure", obj_measure,
                    (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateObjCommand(interp,  "rawtimestep", cmd_rawtimestep,
                    (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateObjCommand(interp,  "gettimestep", cmd_gettimestep,
                    (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

#ifdef VMDPYTHON
  Tcl_CreateCommand(interp, "gopython", text_cmd_gopython,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);
#endif

#if defined(VMDTKCON)
  Tcl_CreateObjCommand(interp,"vmdcon",tcl_vmdcon,
        (ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);
#endif
  
#if defined(VMDCOLVARS)
  Tcl_CreateCommand (interp, "colvars", tcl_colvars, (ClientData) app, (Tcl_CmdDeleteProc*) NULL);
  Tcl_PkgProvide (interp, "colvars", COLVARS_VERSION);
#endif

  Tcl_CreateObjCommand(interp,  "volmap", obj_volmap,
                    (ClientData) app, (Tcl_CmdDeleteProc *) NULL);

  Tcl_CreateCommand(interp, "parallel", text_cmd_parallel,
        (ClientData) app, (Tcl_CmdDeleteProc *) NULL);
  
  return TCL_OK;
}

