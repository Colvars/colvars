// -*- c++ -*-

#include "GlobalMaster.h"
#include "GlobalMasterColvars.h"
#include "colvarproxy_namd.h"


GlobalMasterColvars::GlobalMasterColvars() : proxy(new colvarproxy_namd(this)) {}

GlobalMasterColvars::~GlobalMasterColvars() { reset(); }

void GlobalMasterColvars::calculate() { proxy->calculate(); }


void GlobalMasterColvars::reset()
{
  // Unrequest all positions, total forces, etc from NAMD
  modifyRequestedAtoms().clear();
  modifyForcedAtoms().clear();
  modifyAppliedForces().clear();

  modifyRequestedGroups().clear();
  modifyGroupForces().clear();

#if NAMD_VERSION_NUMBER >= 34471681
  modifyRequestedGridObjects().clear();
  modifyGridObjForces().clear();
#endif

  requestTotalForce(false);
}
