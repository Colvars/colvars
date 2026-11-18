// Adapted from Colvars PR 979. See ../UPSTREAM.md and ../LICENSE.
#ifndef OPENMM_COLVARS_FORCE_H
#define OPENMM_COLVARS_FORCE_H

#include "openmm/Force.h"
#include <string>
#include <vector>

namespace OpenMM { class Context; }
namespace ColvarsPlugin {

/** Colvars biases as an ordinary OpenMM Force.
 * Configuration units are gromacs: nm, kJ/mol, and kJ/mol/nm. Atom numbers
 * in a Colvars configuration are one based. Configure before creating Context.
 */
class ColvarsForce : public OpenMM::Force {
public:
    enum class ConfigSource { File, String };
    enum class PbcMode { Auto, Enabled, Disabled };
    explicit ColvarsForce(const std::string& configFile,
                          const std::string& outputPrefix = "colvars");
    static ColvarsForce fromConfigString(const std::string& config,
                                         const std::string& outputPrefix = "colvars");
    ConfigSource getConfigSource() const { return source_; }
    const std::string& getConfig() const { return config_; }
    const std::string& getOutputPrefix() const { return outputPrefix_; }
    const std::string& getRestartInputPrefix() const { return restartInputPrefix_; }
    void setRestartInputPrefix(const std::string& prefix) { restartInputPrefix_ = prefix; }
    const std::string& getUnitSystem() const;
    void setTargetTemperature(double temperature);
    double getTargetTemperature() const { return temperature_; }
    void setRandomSeed(int seed);
    int getRandomSeed() const { return seed_; }
    void setUseGpu(bool enabled) { useGpu_ = enabled; }
    bool getUseGpu() const { return useGpu_; }
    void setUseGPU(bool enabled) { setUseGpu(enabled); }
    bool getUseGPU() const { return getUseGpu(); }
    /// Fixed charges in elementary-charge units, in zero-based System order.
    /// Required for charge-dependent CVs; supply before constructing Context.
    void setParticleCharges(const std::vector<double>& charges);
    double getParticleCharge(int index) const;
    int getNumParticleCharges() const { return static_cast<int>(charges_.size()); }
    void setPbcMode(PbcMode mode);
    PbcMode getPbcMode() const { return pbcMode_; }
    void setUsesPeriodicBoundaryConditions(bool enabled);
    bool usesPeriodicBoundaryConditions() const override;

    /// Diagnostics describe the last evaluation; these do not advance a bias.
    std::string getBackend(OpenMM::Context& context) const;
    double getColvarValue(OpenMM::Context& context, const std::string& name) const;
    long long getBiasUpdateCount(OpenMM::Context& context) const;
    std::string getBiasState(OpenMM::Context& context) const;

    /// Write the current restart state and flush pending bias output files.
    /// Call after the final energy/force evaluation and before deleting Context.
    /// This does not advance adaptive history, and the simulation may continue.
    void flushOutputs(OpenMM::Context& context) const;

protected:
    OpenMM::ForceImpl* createImpl() const override;
private:
    ConfigSource source_ = ConfigSource::File;
    std::string config_, outputPrefix_, restartInputPrefix_;
    double temperature_ = 300.0;
    int seed_ = 0;
    bool useGpu_ = true;
    PbcMode pbcMode_ = PbcMode::Auto;
    std::vector<double> charges_;
};
}
#endif
