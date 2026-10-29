#include "ColvarsTABFIntegrator.h"

#include "openmm/OpenMMException.h"

#include <algorithm>
#include <cmath>
#include <set>

using namespace OpenMM;

namespace {

void validateForceGroup(int group, const std::string& name)
{
    if (group < 0 || group > 31) {
        throw OpenMMException(name + " must be between 0 and 31");
    }
}

void validateName(const std::string& name, const std::string& label)
{
    if (name.empty()) {
        throw OpenMMException(label + " must not be empty");
    }
}

double validateStepSize(double stepSize)
{
    if (!std::isfinite(stepSize) || stepSize < 0.0) {
        throw OpenMMException("stepSize must be finite and non-negative");
    }
    return stepSize;
}

void validateFiniteNonnegative(double value, const std::string& label)
{
    if (!std::isfinite(value) || value < 0.0) {
        throw OpenMMException(label + " must be finite and non-negative");
    }
}

std::string forceVariable(int group)
{
    return "f" + std::to_string(group);
}

std::string energyVariable(int group)
{
    return "energy" + std::to_string(group);
}

} // namespace

ColvarsTABFIntegrator::ColvarsTABFIntegrator(double stepSize,
                                             double temperature,
                                             double frictionCoeff,
                                             double initialAlpha,
                                             int physicalForceGroup,
                                             int colvarsForceGroup,
                                             const std::string& alphaVariableName,
                                             const std::string& alphaDerivativeVariableName)
    : CustomIntegrator(validateStepSize(stepSize)),
      physicalForceGroup_(physicalForceGroup),
      colvarsForceGroup_(colvarsForceGroup),
      alphaVariableName_(alphaVariableName),
      alphaDerivativeVariableName_(alphaDerivativeVariableName)
{
    validateForceGroup(physicalForceGroup_, "physicalForceGroup");
    validateForceGroup(colvarsForceGroup_, "colvarsForceGroup");
    if (physicalForceGroup_ == colvarsForceGroup_) {
        throw OpenMMException("physicalForceGroup and colvarsForceGroup must be different");
    }
    validateName(alphaVariableName_, "alphaVariableName");
    validateName(alphaDerivativeVariableName_, "alphaDerivativeVariableName");
    if (alphaVariableName_ == alphaDerivativeVariableName_) {
        throw OpenMMException("alpha variable and derivative variable names must be different");
    }
    static const std::set<std::string> reservedNames = {"a", "b", "kT", "x1"};
    if (reservedNames.count(alphaVariableName_) != 0 ||
        reservedNames.count(alphaDerivativeVariableName_) != 0) {
        throw OpenMMException("alpha variable names must not conflict with a, b, kT, or x1");
    }
    validateFiniteNonnegative(temperature, "temperature");
    validateFiniteNonnegative(frictionCoeff, "frictionCoeff");
    if (!std::isfinite(initialAlpha) || initialAlpha <= 0.0) {
        throw OpenMMException("initialAlpha must be finite and positive");
    }

    const double a = std::exp(-frictionCoeff * stepSize);
    const double b = std::sqrt(std::max(0.0, 1.0 - std::exp(-2.0 * frictionCoeff * stepSize)));

    addGlobalVariable("a", a);
    addGlobalVariable("b", b);
    addGlobalVariable("kT", 0.00831446261815324 * temperature);
    addGlobalVariable(alphaVariableName_, initialAlpha);
    addGlobalVariable(alphaDerivativeVariableName_, 0.0);
    addPerDofVariable("x1", 0.0);

    const std::string physicalForce = forceVariable(physicalForceGroup_);
    const std::string colvarsForce = forceVariable(colvarsForceGroup_);
    const std::string physicalEnergy = energyVariable(physicalForceGroup_);

    addUpdateContextState();

    addComputePerDof("v", "v+0.5*dt*" + physicalForce + "/(m*" + alphaVariableName_ + ")");
    addComputePerDof("x", "x+0.5*dt*v");
    addComputePerDof("v", "a*v + b*sqrt(kT/m)*gaussian");
    addComputePerDof("x", "x+0.5*dt*v");
    addComputePerDof("x1", "x");
    addConstrainPositions();
    addComputePerDof("v", "v+(x-x1)/dt");
    addConstrainVelocities();

    addComputeGlobal(alphaDerivativeVariableName_,
                     "-" + physicalEnergy + "/(" + alphaVariableName_ + "*" +
                     alphaVariableName_ + ")");
    addComputePerDof("v", "v+0.5*dt*" + physicalForce + "/(m*" + alphaVariableName_ + ")");
    // Colvars owns history-dependent state (ABF counts, trajectories, alpha).
    // Evaluate it once at the end of the MD step, then apply the full impulse.
    addComputePerDof("v", "v+dt*" + colvarsForce + "/m");
    addConstrainVelocities();
}

int ColvarsTABFIntegrator::getPhysicalForceGroup() const
{
    return physicalForceGroup_;
}

int ColvarsTABFIntegrator::getColvarsForceGroup() const
{
    return colvarsForceGroup_;
}

const std::string& ColvarsTABFIntegrator::getAlphaVariableName() const
{
    return alphaVariableName_;
}

const std::string& ColvarsTABFIntegrator::getAlphaDerivativeVariableName() const
{
    return alphaDerivativeVariableName_;
}
