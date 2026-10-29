#ifndef OPENMM_COLVARS_TABF_INTEGRATOR_H
#define OPENMM_COLVARS_TABF_INTEGRATOR_H

#include "openmm/CustomIntegrator.h"

#include <string>

namespace OpenMM {

/// Langevin-middle style integrator for Colvars TABF simulations.
///
/// This integrator keeps the physical potential scaled as U/alpha by applying
/// physical forces as f/alpha, while applying Colvars bias forces unscaled.
/// Colvars owns the alpha extended-variable dynamics through tabfAlpha; this
/// integrator publishes d(U/alpha)/dalpha = -U/alpha^2 from the physical force
/// group for Colvars to consume.
class ColvarsTABFIntegrator : public CustomIntegrator {
public:
    ColvarsTABFIntegrator(double stepSize,
                          double temperature,
                          double frictionCoeff,
                          double initialAlpha = 1.0,
                          int physicalForceGroup = 0,
                          int colvarsForceGroup = 1,
                          const std::string& alphaVariableName = "tabf_alpha",
                          const std::string& alphaDerivativeVariableName = "tabf_dE_dalpha");

    int getPhysicalForceGroup() const;
    int getColvarsForceGroup() const;
    const std::string& getAlphaVariableName() const;
    const std::string& getAlphaDerivativeVariableName() const;

private:
    int physicalForceGroup_;
    int colvarsForceGroup_;
    std::string alphaVariableName_;
    std::string alphaDerivativeVariableName_;
};

} // namespace OpenMM

#endif
