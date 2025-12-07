#ifndef GLOBALMASTERCOLVARS_H
#define GLOBALMASTERCOLVARS_H

#include <memory>

#include "GlobalMaster.h"


class colvarproxy_namd;

class GlobalMasterColvars : public GlobalMaster {
public:

  GlobalMasterColvars();

  ~GlobalMasterColvars();

  void reset();

  void calculate() override;

  inline void requestTotalForcePublic(bool yesno = true)
  {
    requestTotalForce(yesno);
  }

  inline AtomIDList &modifyRequestedAtomsPublic()
  {
    return modifyRequestedAtoms();
  }

  inline AtomIDList &modifyForcedAtomsPublic()
  {
    return modifyForcedAtoms();
  }

  inline ForceList &modifyAppliedForcesPublic()
  {
    return modifyAppliedForces();
  }

  inline ResizeArray<AtomIDList> const &getRequestedGroups() const
  {
    return reqGroups;
  }

  inline ResizeArray<AtomIDList> &modifyRequestedGroupsPublic()
  {
    return modifyRequestedGroups();
  }

  inline ForceList &modifyGroupForcesPublic()
  {
    return modifyGroupForces();
  }

  inline IntList const &getRequestedGridObjects() const
  {
    return reqGridObjs;
  }

  inline IntList &modifyRequestedGridObjectsPublic()
  {
    return modifyRequestedGridObjects();
  }

  inline BigRealList &modifyGridObjForcesPublic()
  {
    return modifyGridObjForces();
  }

  inline AtomIDList::const_iterator getAtomIdBeginPublic()
  {
    return getAtomIdBegin();
  }

  inline AtomIDList::const_iterator getAtomIdEndPublic()
  {
    return getAtomIdEnd();
  }

  inline PositionList::const_iterator getAtomPositionBeginPublic()
  {
    return getAtomPositionBegin();
  }

  inline PositionList::const_iterator getGroupPositionBeginPublic()
  {
    return getGroupPositionBegin();
  }

  inline PositionList::const_iterator getGroupPositionEndPublic()
  {
    return getGroupPositionEnd();
  }

  inline ForceList::const_iterator getGroupTotalForceBeginPublic()
  {
    return getGroupTotalForceBegin();
  }

  inline ForceList::const_iterator getGroupTotalForceEndPublic()
  {
    return getGroupTotalForceEnd();
  }

  inline IntList::const_iterator getGridObjIndexBeginPublic()
  {
    return getGridObjIndexBegin();
  }

  inline IntList::const_iterator getGridObjIndexEndPublic()
  {
    return getGridObjIndexEnd();
  }

  inline BigRealList::const_iterator getGridObjValueBeginPublic()
  {
    return getGridObjValueBegin();
  }

  inline BigRealList::const_iterator getGridObjValueEndPublic()
  {
    return getGridObjValueEnd();
  }

  inline AtomIDList::const_iterator getForceIdBeginPublic()
  {
    return getForceIdBegin();
  }

  inline AtomIDList::const_iterator getForceIdEndPublic()
  {
    return getForceIdEnd();
  }

  inline ForceList::const_iterator getTotalForcePublic()
  {
    return getTotalForce();
  }

  inline void addReductionEnergyPublic(int reductionTag, BigReal energy)
  {
    addReductionEnergy(reductionTag, energy);
  }

  inline Lattice const *get_lattice() const
  {
    return lattice;
  }

protected:

  std::unique_ptr<colvarproxy_namd> proxy;
};


namespace {
  // Constants used for profiling CkLoop calls
  constexpr int32_t GLOBAL_MASTER_CKLOOP_CALC_ITEM = 2000;
  constexpr int32_t GLOBAL_MASTER_CKLOOP_CALC_BIASES = 2001;
  constexpr int32_t GLOBAL_MASTER_CKLOOP_CALC_SCRIPTED_BIASES = 2002;
}

#endif
