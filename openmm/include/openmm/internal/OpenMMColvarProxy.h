#ifndef OPENMM_INTERNAL_OPENMM_COLVAR_PROXY_H
#define OPENMM_INTERNAL_OPENMM_COLVAR_PROXY_H

#include "openmm/ColvarsForce.h"
#include "openmm/Vec3.h"

#include "colvarproxy.h"

#include <memory>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace OpenMM {

class ContextImpl;
class System;

/// Colvars proxy implementation backed by an OpenMM Context.
///
/// The proxy owns no OpenMM objects. Context and position pointers are valid
/// only during the current computeForce() call.
class OpenMMColvarProxy final : public colvarproxy {
public:
    OpenMMColvarProxy(const ColvarsForce& owner, const System& system);
    ~OpenMMColvarProxy() override;

    int setup() override;

    int init_atom(int atomNumber) override;
    void clear_atom(int index) override;
    int check_atom_id(int atomNumber) override;

    void request_total_force(bool yesno) override;
    bool total_forces_enabled() const override;
    bool total_forces_same_step() const override;

    int set_unit_system(std::string const& unitsIn, bool checkOnly) override;

    void add_energy(cvm::real energy) override;
    void log(std::string const& message) override;
    void error(std::string const& message) override;
    cvm::real rand_gaussian() override;

    int get_alch_lambda(cvm::real* lambda) override;
    int send_alch_lambda() override;
    int request_alch_energy_freq(int freq) override;
    int get_dE_dlambda(cvm::real* dEdLambda) override;

#ifdef OPENMM_COLVARS_HAS_TABF_ALPHA
    int get_tabf_alpha(cvm::real* alpha) override;
#ifdef OPENMM_COLVARS_HAS_TABF_ALPHA_BOUNDARIES
    int get_tabf_alpha_boundaries(cvm::real* lower, cvm::real* upper) override;
#endif
    int send_tabf_alpha() override;
    int request_tabf_energy_freq(int freq) override;
    int get_dE_dtabf_alpha(cvm::real* dEdAlpha) override;
#endif

    int load_atoms_pdb(char const* filename,
                       cvm::atom_group& atoms,
                       std::string const& pdbField,
                       double pdbFieldValue) override;
    int load_coords_pdb(char const* filename,
                        std::vector<cvm::atom_pos>& positions,
                        std::vector<int> const& sortedIds,
                        std::string const& pdbField,
                        double pdbFieldValue) override;

    int check_replicas_enabled() override;
    int replica_index() override;
    int num_replicas() override;
    void replica_comm_barrier() override;
    int replica_comm_recv(char* msgData, int bufferLength, int sourceReplica) override;
    int replica_comm_send(char* msgData, int messageLength, int destinationReplica) override;

    int run_force_callback() override;
    int run_colvar_callback(std::string const& name,
                            std::vector<const colvarvalue*> const& cvcs,
                            colvarvalue& value) override;
    int run_colvar_gradient_callback(std::string const& name,
                                     std::vector<const colvarvalue*> const& cvcs,
                                     std::vector<cvm::matrix2d<cvm::real>>& gradient) override;

    int updateInput(ContextImpl& context, const std::vector<Vec3>& positions);
    double calculate(ContextImpl& context, const std::vector<Vec3>& positions);
    void addForcesToOpenMM(std::vector<Vec3>& forces) const;

    double biasEnergy() const;

private:
    void initializeColvars();
    void updateAtomProperties(int proxyIndex);
    void updateAtomIndexMap();
    void updatePbc(ContextImpl& context);
    bool shouldUsePbc(ContextImpl& context) const;
    bool resolveAutoPbcMode() const;
    bool pbcBoxChanged(const Vec3& a, const Vec3& b, const Vec3& c) const;
    bool buildCacheKey(ContextImpl& context,
                       const std::vector<Vec3>& positions,
                       long long& openmmStep,
                       bool& usesPbc,
                       Vec3& boxA,
                       Vec3& boxB,
                       Vec3& boxC,
                       double& lambdaValue,
                       double& tabfAlphaValue,
                       std::vector<Vec3>& requestedPositions) const;
    bool cacheMatches(long long openmmStep,
                      bool usesPbc,
                      const Vec3& boxA,
                      const Vec3& boxB,
                      const Vec3& boxC,
                      double lambdaValue,
                      double tabfAlphaValue,
                      const std::vector<Vec3>& requestedPositions) const;
    void storeCache(long long openmmStep,
                    bool usesPbc,
                    const Vec3& boxA,
                    const Vec3& boxB,
                    const Vec3& boxC,
                    double lambdaValue,
                    double tabfAlphaValue);
    void invalidateCache();
    void updateColvarsStep(ContextImpl& context);
    void updatePositions(const std::vector<Vec3>& positions);
    void updateTotalForces(ContextImpl& context);
    bool inferTotalForcesSameStep() const;
    void clearPreviousOpenMMForces(std::vector<Vec3>& forces) const;
    void captureLambdaDerivative(ContextImpl& context);
    double contextParameter(ContextImpl& context, const std::string& name, double defaultValue) const;
    void setContextParameter(ContextImpl& context, const std::string& name, double value) const;
    double integratorGlobal(ContextImpl& context, const std::string& name, double defaultValue) const;
    bool setIntegratorGlobal(ContextImpl& context, const std::string& name, double value) const;
    double tabfAlpha(ContextImpl& context) const;

    const ColvarsForce& owner_;
    const System& system_;
    ContextImpl* currentContext_;
    double biasEnergy_;
    bool isInitialized_;
    mutable bool autoPbcResolved_;
    mutable bool cachedAutoPbc_;
    bool pbcLatticeInitialized_;
    bool reportedRepeatedExtendedCoordinate_;
    bool totalForcesSameStep_;
    bool stepMappingInitialized_;
    long long openmmInitialStep_;
    cvm::step_number colvarsInitialStep_;
    long long lastLambdaWritebackStep_;
    long long lastTabfAlphaWritebackStep_;
    Vec3 lastPbcA_;
    Vec3 lastPbcB_;
    Vec3 lastPbcC_;
    double lastIntegrationTimestepFs_;
    std::unordered_map<int, int> atomIndexById_;
    mutable std::vector<int> previousForceAtoms_;
    std::vector<Vec3> systemForceBuffer_;
    std::uint64_t gaussianRngState_;
    bool gaussianRngHasSpare_;
    cvm::real gaussianRngSpare_;
    bool currentLambdaDerivativeValid_;
    double currentLambdaDerivative_;
    bool tabfAlphaBoundariesValid_;
    cvm::real tabfAlphaLowerBoundary_;
    cvm::real tabfAlphaUpperBoundary_;

    bool cacheValid_;
    long long cachedOpenmmStep_;
    bool cachedUsesPbc_;
    Vec3 cachedPbcA_;
    Vec3 cachedPbcB_;
    Vec3 cachedPbcC_;
    double cachedLambdaValue_;
    double cachedTabfAlphaValue_;
    double cachedBiasEnergy_;
    std::vector<int> cachedAtomIds_;
    std::vector<Vec3> cachedRequestedPositions_;
    std::vector<cvm::rvector> cachedAppliedForces_;
};

} // namespace OpenMM

#endif
