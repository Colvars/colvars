/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/* A ComputeGlobalMaster represents a bit of computation that needs to
   be done on atoms or groups of atoms on several nodes.  It is given
   the positions of atoms and groups, and provides a list of requested
   atoms, forces, and groups in return.
   
   I'm not going to do groups for now, because they were done badly
   originally.  A better solution is necessary.  Hint: multiple
   groups are required.

   The expected usage (by ComputeGlobalMasterServer) is:
   1) receive and store the data message from all nodes
   2) call processData with that data
   3) collect new forces and atoms
   4) send them as a message to all the nodes
   Repeat until done.
*/

#ifndef GLOBALMASTER_H
#define GLOBALMASTER_H

#include "NamdTypes.h"
class Lattice;

class GlobalMaster {
 public:

  /* This passes the atom and group data to the master, which then
     performs any necessary calculations and updates its results. */
  void processData(AtomIDList::iterator a_i,
		   AtomIDList::iterator a_e,
		   PositionList::iterator p_i,
		   PositionList::iterator g_i,
		   PositionList::iterator g_e,
		   BigRealList::iterator gm_i,
		   BigRealList::iterator gm_e,
		   ForceList::iterator gtf_i,
		   ForceList::iterator gtf_e,
                   IntList::iterator goi_i,
                   IntList::iterator goi_e,
                   BigRealList::iterator gov_i,
                   BigRealList::iterator gov_e,
		   AtomIDList::iterator last_atoms_forced_i,
		   AtomIDList::iterator last_atoms_forced_e,
		   ForceList::iterator last_forces_i,
		   AtomIDList::iterator,
		   AtomIDList::iterator,
		   ForceList::iterator);

  int step;  // set by server to current timestep before processData
  int globalMasterStep;  // set by server to current timestep / globalMasterFreq before processData
  int old_num_groups_requested;  // used for group forces

  bool changedAtoms(); // false if the atom IDs haven't changed
  const AtomIDList &requestedAtoms(); // the atom ids requested
  bool changedForces(); // false if the forces haven't changed
  const AtomIDList &forcedAtoms(); // the atoms the forces are applied to
  const ForceList &appliedForces(); // the corresponding applied forces
  bool changedGroups(); // false if the groups haven't changed
  const ResizeArray<AtomIDList> &requestedGroups(); // the requested groups
  const ForceList &groupForces(); // the corresponding forces on groups
  bool changedGridObjs(); // false if the groups haven't changed
  const IntList &requestedGridObjs(); // the requested groups
  const BigRealList &gridObjForces(); // the corresponding forces on groups
  bool requestedTotalForces() { return totalForceRequested; }

  /* sets changedAtoms and changedForces to false again */
  void clearChanged(); 
  virtual ~GlobalMaster() {}; // necessary for abstract classes '-P

  void check() const; // dies if there are problems with the rep invariants

  void setLattice(const Lattice *lat) { lattice = lat; }
  
 protected:
  GlobalMaster();

  /* This will be called after the pointers to lists below have been
     initialized correctly by processData.  It should perform any
     required caluation and update the atom/force lists. */
  virtual void calculate();

  /* This function returns the list of requested atoms, but assumes
     that you will change it. */
  AtomIDList &modifyRequestedAtoms();

  /* These functions returns the list of requested forces, but assumes
     that you will change it.  The two lists must be kept at the same
     length, since the forcedAtoms correspond directly to the
     appliedForces. */
  AtomIDList &modifyForcedAtoms();
  ForceList &modifyAppliedForces();

  /* This function lets you change the requested groups */
  ResizeArray<AtomIDList> &modifyRequestedGroups();
  ForceList &modifyGroupForces();

  /* Same here for grids */
  IntList &modifyRequestedGridObjects();
  BigRealList &modifyGridObjForces();

  /* These return pointers to the lists of atom ids and positions, as
     they were last passed to processData (see below) */
  AtomIDList::const_iterator getAtomIdBegin();
  AtomIDList::const_iterator getAtomIdEnd();
  PositionList::const_iterator getAtomPositionBegin();
  PositionList::const_iterator getGroupPositionBegin();
  PositionList::const_iterator getGroupPositionEnd();
  ForceList::const_iterator getGroupTotalForceBegin();
  ForceList::const_iterator getGroupTotalForceEnd();
  IntList::const_iterator getGridObjIndexBegin();
  IntList::const_iterator getGridObjIndexEnd();
  BigRealList::const_iterator getGridObjValueBegin();
  BigRealList::const_iterator getGridObjValueEnd();
  
  /* these give you all the global forces being applied by masters */
  /* again, here we only need one end iterator */
  AtomIDList::const_iterator getLastAtomsForcedBegin();
  AtomIDList::const_iterator getLastAtomsForcedEnd();
  ForceList::const_iterator getLastForcesBegin();
  
  /* These return the pointers to the lists of requested atom IDs
     and total forces on these atoms */
  AtomIDList::const_iterator getForceIdBegin();
  AtomIDList::const_iterator getForceIdEnd();
  ForceList::const_iterator getTotalForce();

  bool totalForceRequested;
  void requestTotalForce(bool yesno = true) { totalForceRequested = yesno; }
  
  /* This helpful function returns an array with the masses of each of
     the groups whose positions we have.  */
  BigRealList::const_iterator getGroupMassBegin();
  BigRealList::const_iterator getGroupMassEnd();

 protected:
  const Lattice *lattice;  // points to lattice in server

  /* These store the pointers to lists of atom ids and atom positions.
     The list of atom positions has the same length as the list of
     ids, so only three iterators are necessary.   There are also
     pointers to the beginning and end of the group position list
     here. */
  AtomIDList::iterator atomIdBegin;
  AtomIDList::iterator atomIdEnd;
  PositionList::iterator atomPositionBegin;
  PositionList::iterator groupPositionBegin;
  PositionList::iterator groupPositionEnd;
  BigRealList::iterator groupMassBegin;
  BigRealList::iterator groupMassEnd;
  ForceList::iterator groupTotalForceBegin;
  ForceList::iterator groupTotalForceEnd;
  IntList::iterator gridObjIndexBegin;
  IntList::iterator gridObjIndexEnd;
  BigRealList::iterator gridObjValueBegin;
  BigRealList::iterator gridObjValueEnd;

  /* these store all the global forces being applied by masters */
  AtomIDList::iterator lastAtomsForcedBegin;
  ForceList::iterator lastForcesBegin;
  AtomIDList::iterator lastAtomsForcedEnd;
  
  /* These store all the total forces returned from the simulation */
  AtomIDList::iterator forceIdBegin;
  AtomIDList::iterator forceIdEnd;
  ForceList::iterator totalForceBegin;

  /* These store the requested atoms and forces, and the booleans
     indicate whether they (may) have changed. */
  bool reqAtomsChanged;
  AtomIDList reqAtoms; // atoms whose positions are requested

  bool appForcesChanged;
  AtomIDList fAtoms; // atoms that are being forced
  ForceList appForces; // the corresponding forces

  bool reqGroupsChanged;
  ResizeArray<AtomIDList> reqGroups; // list of requested groups of atoms 
  ForceList grpForces; // the corresponding forces

  bool reqGridObjsChanged;
  IntList reqGridObjs; // list of requested grids
  BigRealList gridobjForces; // the corresponding forces

  friend class colvarproxy_namd;
};

#endif
