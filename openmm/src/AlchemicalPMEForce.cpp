#include "AlchemicalPMEForce.h"
#include "AlchemicalPMEForceImpl.h"

#include "openmm/OpenMMException.h"

#include <cmath>

using namespace OpenMM;

namespace {

void validateIndex(int index, int size, const std::string& name)
{
    if (index < 0 || index >= size) {
        throw OpenMMException(name + " index out of range");
    }
}

void validateName(const std::string& name, const std::string& label)
{
    if (name.empty()) {
        throw OpenMMException(label + " must not be empty");
    }
}

void validateNonnegative(double value, const std::string& label)
{
    if (!std::isfinite(value) || value < 0.0) {
        throw OpenMMException(label + " must be finite and non-negative");
    }
}

void validateFinite(double value, const std::string& label)
{
    if (!std::isfinite(value)) {
        throw OpenMMException(label + " must be finite");
    }
}

} // namespace

AlchemicalPMEForce::AlchemicalPMEForce()
    : lambdaParameterName_("lambda"),
      lambdaDerivativeParameterName_("dE_dlambda"),
      electrostaticLambdaStart_(0.0),
      nonbondedMethod_(PME),
      cutoffDistance_(1.0),
      reactionFieldDielectric_(78.3),
      pmeAlpha_(0.0),
      pmeGridX_(0),
      pmeGridY_(0),
      pmeGridZ_(0),
      pmeInterpolationOrder_(5)
{
}

int AlchemicalPMEForce::addParticle(double charge0, double charge1,
                                    double sigma0, double sigma1,
                                    double epsilon0, double epsilon1)
{
    validateFinite(charge0, "charge0");
    validateFinite(charge1, "charge1");
    validateNonnegative(sigma0, "sigma0");
    validateNonnegative(sigma1, "sigma1");
    validateNonnegative(epsilon0, "epsilon0");
    validateNonnegative(epsilon1, "epsilon1");
    particles_.push_back({charge0, charge1, sigma0, sigma1, epsilon0, epsilon1});
    return static_cast<int>(particles_.size()) - 1;
}

int AlchemicalPMEForce::getNumParticles() const
{
    return static_cast<int>(particles_.size());
}

void AlchemicalPMEForce::getParticleParameters(int index,
                                               double& charge0, double& charge1,
                                               double& sigma0, double& sigma1,
                                               double& epsilon0, double& epsilon1) const
{
    validateIndex(index, getNumParticles(), "particle");
    const auto& p = particles_[index];
    charge0 = p.charge0;
    charge1 = p.charge1;
    sigma0 = p.sigma0;
    sigma1 = p.sigma1;
    epsilon0 = p.epsilon0;
    epsilon1 = p.epsilon1;
}

void AlchemicalPMEForce::setParticleParameters(int index,
                                               double charge0, double charge1,
                                               double sigma0, double sigma1,
                                               double epsilon0, double epsilon1)
{
    validateIndex(index, getNumParticles(), "particle");
    validateFinite(charge0, "charge0");
    validateFinite(charge1, "charge1");
    validateNonnegative(sigma0, "sigma0");
    validateNonnegative(sigma1, "sigma1");
    validateNonnegative(epsilon0, "epsilon0");
    validateNonnegative(epsilon1, "epsilon1");
    particles_[index] = {charge0, charge1, sigma0, sigma1, epsilon0, epsilon1};
}

int AlchemicalPMEForce::addException(int particle1, int particle2,
                                     double chargeProd0, double chargeProd1,
                                     double sigma0, double sigma1,
                                     double epsilon0, double epsilon1)
{
    validateIndex(particle1, getNumParticles(), "exception particle1");
    validateIndex(particle2, getNumParticles(), "exception particle2");
    if (particle1 == particle2) {
        throw OpenMMException("An AlchemicalPMEForce exception cannot connect a particle to itself");
    }
    validateFinite(chargeProd0, "chargeProd0");
    validateFinite(chargeProd1, "chargeProd1");
    validateNonnegative(sigma0, "sigma0");
    validateNonnegative(sigma1, "sigma1");
    validateNonnegative(epsilon0, "epsilon0");
    validateNonnegative(epsilon1, "epsilon1");
    exceptions_.push_back({particle1, particle2, chargeProd0, chargeProd1, sigma0, sigma1, epsilon0, epsilon1});
    return static_cast<int>(exceptions_.size()) - 1;
}

int AlchemicalPMEForce::getNumExceptions() const
{
    return static_cast<int>(exceptions_.size());
}

void AlchemicalPMEForce::getExceptionParameters(int index,
                                                int& particle1, int& particle2,
                                                double& chargeProd0, double& chargeProd1,
                                                double& sigma0, double& sigma1,
                                                double& epsilon0, double& epsilon1) const
{
    validateIndex(index, getNumExceptions(), "exception");
    const auto& e = exceptions_[index];
    particle1 = e.particle1;
    particle2 = e.particle2;
    chargeProd0 = e.chargeProd0;
    chargeProd1 = e.chargeProd1;
    sigma0 = e.sigma0;
    sigma1 = e.sigma1;
    epsilon0 = e.epsilon0;
    epsilon1 = e.epsilon1;
}

void AlchemicalPMEForce::setExceptionParameters(int index,
                                                int particle1, int particle2,
                                                double chargeProd0, double chargeProd1,
                                                double sigma0, double sigma1,
                                                double epsilon0, double epsilon1)
{
    validateIndex(index, getNumExceptions(), "exception");
    validateIndex(particle1, getNumParticles(), "exception particle1");
    validateIndex(particle2, getNumParticles(), "exception particle2");
    if (particle1 == particle2) {
        throw OpenMMException("An AlchemicalPMEForce exception cannot connect a particle to itself");
    }
    validateFinite(chargeProd0, "chargeProd0");
    validateFinite(chargeProd1, "chargeProd1");
    validateNonnegative(sigma0, "sigma0");
    validateNonnegative(sigma1, "sigma1");
    validateNonnegative(epsilon0, "epsilon0");
    validateNonnegative(epsilon1, "epsilon1");
    exceptions_[index] = {particle1, particle2, chargeProd0, chargeProd1, sigma0, sigma1, epsilon0, epsilon1};
}

const std::string& AlchemicalPMEForce::getLambdaParameterName() const
{
    return lambdaParameterName_;
}

void AlchemicalPMEForce::setLambdaParameterName(const std::string& name)
{
    validateName(name, "lambdaParameterName");
    lambdaParameterName_ = name;
}

const std::string& AlchemicalPMEForce::getLambdaDerivativeParameterName() const
{
    return lambdaDerivativeParameterName_;
}

void AlchemicalPMEForce::setLambdaDerivativeParameterName(const std::string& name)
{
    validateName(name, "lambdaDerivativeParameterName");
    lambdaDerivativeParameterName_ = name;
}

double AlchemicalPMEForce::getElectrostaticLambdaStart() const
{
    return electrostaticLambdaStart_;
}

void AlchemicalPMEForce::setElectrostaticLambdaStart(double lambdaStart)
{
    if (!std::isfinite(lambdaStart) || lambdaStart < 0.0 || lambdaStart >= 1.0) {
        throw OpenMMException("electrostaticLambdaStart must be finite and in the range [0, 1)");
    }
    electrostaticLambdaStart_ = lambdaStart;
}

AlchemicalPMEForce::NonbondedMethod AlchemicalPMEForce::getNonbondedMethod() const
{
    return nonbondedMethod_;
}

void AlchemicalPMEForce::setNonbondedMethod(NonbondedMethod method)
{
    switch (method) {
        case NoCutoff:
        case CutoffPeriodic:
        case PME:
            nonbondedMethod_ = method;
            return;
    }
    throw OpenMMException("Invalid AlchemicalPMEForce nonbonded method");
}

double AlchemicalPMEForce::getCutoffDistance() const
{
    return cutoffDistance_;
}

void AlchemicalPMEForce::setCutoffDistance(double distance)
{
    if (!std::isfinite(distance) || distance <= 0.0) {
        throw OpenMMException("cutoffDistance must be finite and positive");
    }
    cutoffDistance_ = distance;
}

double AlchemicalPMEForce::getReactionFieldDielectric() const
{
    return reactionFieldDielectric_;
}

void AlchemicalPMEForce::setReactionFieldDielectric(double dielectric)
{
    if (!std::isfinite(dielectric) || dielectric < 1.0) {
        throw OpenMMException("reactionFieldDielectric must be finite and at least 1");
    }
    reactionFieldDielectric_ = dielectric;
}

void AlchemicalPMEForce::getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const
{
    alpha = pmeAlpha_;
    nx = pmeGridX_;
    ny = pmeGridY_;
    nz = pmeGridZ_;
}

void AlchemicalPMEForce::setPMEParameters(double alpha, int nx, int ny, int nz)
{
    validateNonnegative(alpha, "alpha");
    const bool automatic = alpha == 0.0 && nx == 0 && ny == 0 && nz == 0;
    const bool explicitParameters = alpha > 0.0 && nx > 0 && ny > 0 && nz > 0;
    if (!automatic && !explicitParameters) {
        throw OpenMMException(
            "PME parameters must either all be zero or have positive alpha and grid dimensions");
    }
    pmeAlpha_ = alpha;
    pmeGridX_ = nx;
    pmeGridY_ = ny;
    pmeGridZ_ = nz;
}

int AlchemicalPMEForce::getPMEInterpolationOrder() const
{
    return pmeInterpolationOrder_;
}

void AlchemicalPMEForce::setPMEInterpolationOrder(int order)
{
    if (order < 3 || order > 6) {
        throw OpenMMException("PME interpolation order must be between 3 and 6");
    }
    pmeInterpolationOrder_ = order;
}

bool AlchemicalPMEForce::usesPeriodicBoundaryConditions() const
{
    return nonbondedMethod_ == CutoffPeriodic || nonbondedMethod_ == PME;
}

ForceImpl* AlchemicalPMEForce::createImpl() const
{
    return new AlchemicalPMEForceImpl(*this);
}
