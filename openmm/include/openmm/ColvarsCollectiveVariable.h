#ifndef OPENMM_COLVARS_COLLECTIVE_VARIABLE_H
#define OPENMM_COLVARS_COLLECTIVE_VARIABLE_H

#include "openmm/Force.h"

#include <string>

namespace OpenMM {

/// OpenMM Force whose energy is the value of one Colvars collective variable.
///
/// The energy is the value of the colvar named `colvarName` in the Colvars
/// configuration, and the forces are minus its gradient with respect to the
/// particle positions. The value is in Colvars units (nanometers for lengths,
/// degrees for angles); a customFunction colvar can convert it, for example to
/// radians. The value of a periodic colvar is reported wrapped into its period,
/// so this force is meant to be used as a collective variable of another force,
/// such as a CustomCVForce, which handles the periodicity.
///
/// The configuration may define other colvars, for example to combine them in a
/// customFunction, but no biases, and the reported colvar must be a scalar without
/// an extended Lagrangian. No Colvars output files are written, and logLevel
/// defaults to 0.
///
/// Unlike ColvarsForce, this force can be serialized, so it can be copied and used
/// inside a CustomCVForce.
class ColvarsCollectiveVariable : public Force {
public:
    ColvarsCollectiveVariable(const std::string& configString, const std::string& colvarName);

    const std::string& getConfigString() const;
    const std::string& getColvarName() const;

    bool usesPeriodicBoundaryConditions() const override;
    void setUsesPeriodicBoundaryConditions(bool usesPbc);

protected:
    ForceImpl* createImpl() const override;

private:
    std::string configString_;
    std::string colvarName_;
    bool usesPbc_;
};

} // namespace OpenMM

#endif
