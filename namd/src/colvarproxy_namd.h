/// -*- c++ -*-

#ifndef COLVARPROXY_NAMD_H
#define COLVARPROXY_NAMD_H

#include "Vector.h"
#include "ResizeArray.h"
#include "NamdTypes.h"
#include "SimParameters.h"
#include "Lattice.h"
#include "GlobalMaster.h"
#include "Random.h"

#include "colvarmodule.h"
#include "colvarproxy.h"
#include "colvarvalue.h"

#ifndef COLVARPROXY_VERSION
#define COLVARPROXY_VERSION "2015-07-24"
#endif

// For replica exchange
#include "converse.h"
#include "DataExchanger.h"

/// \brief Communication between colvars and NAMD (implementation of
/// \link colvarproxy \endlink)
class colvarproxy_namd : public colvarproxy, public GlobalMaster {

protected:

  /// Pointer to the NAMD simulation input object
  SimParameters const *simparams;

  BigReal thermostat_temperature;

  /// NAMD-style PRNG object
  Random random;


  size_t      restart_frequency_s;
  size_t      previous_NAMD_step;
  bool        first_timestep;
  bool        system_force_requested;

  std::vector<int>          colvars_atoms;
  std::vector<size_t>       colvars_atoms_ncopies;
  std::vector<cvm::rvector> positions;
  std::vector<cvm::rvector> total_forces;
  std::vector<cvm::rvector> applied_forces;

  size_t init_namd_atom(AtomID const &aid);

  SubmitReduction *reduction;

#ifdef NAMD_TCL
  Tcl_Interp *interp; // Tcl interpreter embedded in NAMD
#endif

public:

  friend class cvm::atom;

  colvarproxy_namd();
  ~colvarproxy_namd();

  /// \brief Reimplements GlobalMaster member function, to be called
  /// at every step
  void calculate();

  void add_energy(cvm::real energy);
  void request_system_force(bool yesno);
  void log(std::string const &message);
  void error(std::string const &message);
  void fatal_error(std::string const &message);
  void exit(std::string const &message);

  // Callback functions
  int run_force_callback();
  int run_colvar_callback(std::string const &name,
                          std::vector<const colvarvalue *> const &cvcs,
                          colvarvalue &value);
  int run_colvar_gradient_callback(std::string const &name,
                                   std::vector<const colvarvalue *> const &cvcs,
                                   std::vector<cvm::matrix2d<cvm::real> > &gradient);

  inline cvm::real unit_angstrom()
  {
    return 1.0;
  }

  cvm::real boltzmann()
  {
    return 0.001987191;
  }

  cvm::real temperature()
  {
    return thermostat_temperature;
  }

  cvm::real dt()
  {
    return simparams->dt;
  }

  // Replica communication functions.
  bool replica_enabled() {
#if CMK_HAS_PARTITION
    return true;
#else
    return false;
#endif
  }

  int replica_index() {
    return CmiMyPartition();
  }

  int replica_num() {
    return CmiNumPartitions();
  }

  void replica_comm_barrier() {
    replica_barrier();
  }

  int replica_comm_recv(char* msg_data, int buf_len, int src_rep) {
    DataMessage *recvMsg = NULL;
    replica_recv(&recvMsg, src_rep, CkMyPe());
    CmiAssert(recvMsg != NULL);
    int retval = recvMsg->size;
    if (buf_len >= retval) {
      memcpy(msg_data,recvMsg->data,retval);
    } else {
      retval = 0;
    }
    CmiFree(recvMsg);
    return retval;
  }

  int replica_comm_send(char* msg_data, int msg_len, int dest_rep) {
    replica_send(msg_data, msg_len, dest_rep, CkMyPe());
    return msg_len;
  }

  inline size_t restart_frequency()
  {
    return restart_frequency_s;
  }

  cvm::rvector position_distance(cvm::atom_pos const &pos1,
                                 cvm::atom_pos const &pos2);
  cvm::real position_dist2(cvm::atom_pos const &pos1,
                           cvm::atom_pos const &pos2);

  void select_closest_image(cvm::atom_pos &pos,
                            cvm::atom_pos const &ref_pos);


  int load_atoms(char const *filename,
                 std::vector<cvm::atom> &atoms,
                 std::string const &pdb_field,
                 double const pdb_field_value = 0.0);

  int load_coords(char const *filename,
                  std::vector<cvm::atom_pos> &pos,
                  const std::vector<int> &indices,
                  std::string const &pdb_field,
                  double const pdb_field_value = 0.0);

  std::ostream * output_stream(std::string const &output_name);
  int close_output_stream(std::string const &output_name);
  int backup_file(char const *filename);

  cvm::real rand_gaussian(void)
  {
    return random.gaussian();
  }
};


inline cvm::rvector colvarproxy_namd::position_distance(cvm::atom_pos const &pos1,
                                                        cvm::atom_pos const &pos2)
{
  Position const p1(pos1.x, pos1.y, pos1.z);
  Position const p2(pos2.x, pos2.y, pos2.z);
  // return p2 - p1
  Vector const d = this->lattice->delta(p2, p1);
  return cvm::rvector(d.x, d.y, d.z);
}


inline void colvarproxy_namd::select_closest_image(cvm::atom_pos &pos,
                                                   cvm::atom_pos const &ref_pos)
{
  Position const p(pos.x, pos.y, pos.z);
  Position const rp(ref_pos.x, ref_pos.y, ref_pos.z);
  ScaledPosition const srp = this->lattice->scale(rp);
  Position const np = this->lattice->nearest(p, srp);
  pos.x = np.x;
  pos.y = np.y;
  pos.z = np.z;
}


inline cvm::real colvarproxy_namd::position_dist2(cvm::atom_pos const &pos1,
                                                  cvm::atom_pos const &pos2)
{
  Lattice const *l = this->lattice;
  Vector const p1(pos1.x, pos1.y, pos1.z);
  Vector const p2(pos2.x, pos2.y, pos2.z);
  Vector const d = l->delta(p1, p2);
  return cvm::real(d.x*d.x + d.y*d.y + d.z*d.z);
}



#endif
