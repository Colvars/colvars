#ifndef COLVARPROXY_CUDAGLOBALMASTER_H
#define COLVARPROXY_CUDAGLOBALMASTER_H

#include "CudaGlobalMasterClient.h"
#include <string>
#include <vector>
#include <memory>
#include <cuda_runtime.h>

class Lattice;

class colvarproxy_impl;

// NOTE: Don't inherit from both CudaGlobalMasterClient and colvarproxy!
class CudaGlobalMasterColvars: public CudaGlobalMasterClient {
public:
  CudaGlobalMasterColvars();
  virtual ~CudaGlobalMasterColvars();
  void initialize(const std::vector<std::string>& arguments, int deviceID, cudaStream_t stream) override;
  void calculate() override;
  cudaStream_t getStream() override;
  bool requestedAtomsChanged() override;
  bool requestedForcedAtomsChanged() override { return requestedAtomsChanged(); }
  bool requestedTotalForcesAtomsChanged() override;
  bool requestUpdateAtomPositions() override { return true; }
  bool requestUpdateAtomTotalForces() override;
  bool requestUpdateForcedAtoms() override { return requestUpdateAtomPositions(); }
  bool requestUpdateMasses() override;
  bool requestUpdateCharges() override;
  bool requestUpdateLattice() override { return true; }
  double getEnergy() const override;
  double* getAppliedForces() const override;
  double* getPositions() override;
  float* getMasses() override;
  float* getCharges() override;
  double* getTotalForces() override;
  double* getLattice() override;
  const std::vector<AtomID>& getRequestedAtoms() const override;
  const std::vector<AtomID>& getRequestedTotalForcesAtoms() const override;
  const std::vector<AtomID>& getRequestedForcedAtoms() const override {
    return getRequestedAtoms();
  }
  bool replica_enabled() const {
    return CudaGlobalMasterClient::replica_enabled();
  }
  int replica_index() const {
    return CudaGlobalMasterClient::replica_index();
  }
  int num_replicas() const {
    return CudaGlobalMasterClient::num_replicas();
  }
  void replica_comm_barrier() {
    CudaGlobalMasterClient::replica_comm_barrier();
  }
  int replica_comm_recv(char* msg_data, int buf_len, int src_rep) {
    return CudaGlobalMasterClient::replica_comm_recv(msg_data, buf_len, src_rep);
  }
  int replica_comm_send(char* msg_data, int msg_len, int dest_rep) {
    return CudaGlobalMasterClient::replica_comm_send(msg_data, msg_len, dest_rep);
  }
  int64_t getStep() const {return m_step;}
private:
  std::unique_ptr<colvarproxy_impl> mImpl;
  std::vector<AtomID> mEmpty;
};

#endif // COLVARPROXY_CUDAGLOBALMASTER_H
