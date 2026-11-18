// Adapted from Colvars PR 979. See ../UPSTREAM.md and ../LICENSE.
#include "ColvarsForce.h"
#include "ColvarsForceImpl.h"
#include "colvarproxy_openmm.h"
#include "openmm/OpenMMException.h"
#include <cmath>
using namespace OpenMM;
using namespace ColvarsPlugin;

ColvarsForce::ColvarsForce(const std::string& configFile, const std::string& outputPrefix)
    : config_(configFile), outputPrefix_(outputPrefix) {
    if (config_.empty() || outputPrefix_.empty())
        throw OpenMMException("Colvars configuration and output prefix must be nonempty");
}
ColvarsForce ColvarsForce::fromConfigString(const std::string& config, const std::string& prefix) {
    ColvarsForce force(config, prefix);
    force.source_ = ConfigSource::String;
    return force;
}
const std::string& ColvarsForce::getUnitSystem() const {
    static const std::string units("gromacs");
    return units;
}
void ColvarsForce::setTargetTemperature(double temperature) {
    if (!std::isfinite(temperature) || temperature <= 0)
        throw OpenMMException("Colvars temperature must be positive and finite");
    temperature_ = temperature;
}
void ColvarsForce::setRandomSeed(int seed) {
    if (seed < 0) throw OpenMMException("Colvars random seed must be nonnegative");
    seed_ = seed;
}
void ColvarsForce::setParticleCharges(const std::vector<double>& charges) {
    for (double q : charges)
        if (!std::isfinite(q)) throw OpenMMException("Colvars particle charges must be finite");
    charges_ = charges;
}
double ColvarsForce::getParticleCharge(int i) const {
    return i >= 0 && i < static_cast<int>(charges_.size()) ? charges_[i] : 0.0;
}
void ColvarsForce::setPbcMode(PbcMode mode) {
    if (mode != PbcMode::Auto && mode != PbcMode::Enabled && mode != PbcMode::Disabled)
        throw OpenMMException("Invalid Colvars PBC mode");
    pbcMode_ = mode;
}
void ColvarsForce::setUsesPeriodicBoundaryConditions(bool enabled) {
    pbcMode_ = enabled ? PbcMode::Enabled : PbcMode::Disabled;
}
bool ColvarsForce::usesPeriodicBoundaryConditions() const { return pbcMode_ == PbcMode::Enabled; }
ForceImpl* ColvarsForce::createImpl() const { return new ColvarsForceImpl(*this); }
std::string ColvarsForce::getBackend(Context& c) const {
    return dynamic_cast<ColvarsForceImpl&>(getImplInContext(c)).backend();
}
double ColvarsForce::getColvarValue(Context& c, const std::string& name) const {
    return dynamic_cast<ColvarsForceImpl&>(getImplInContext(c)).proxy().getColvarValue(name);
}
long long ColvarsForce::getBiasUpdateCount(Context& c) const {
    return dynamic_cast<ColvarsForceImpl&>(getImplInContext(c)).proxy().getBiasUpdateCount();
}
std::string ColvarsForce::getBiasState(Context& c) const {
    return dynamic_cast<ColvarsForceImpl&>(getImplInContext(c)).proxy().getBiasState();
}
void ColvarsForce::flushOutputs(Context& c) const {
    dynamic_cast<ColvarsForceImpl&>(getImplInContext(c)).proxy().flushOutputs();
}
