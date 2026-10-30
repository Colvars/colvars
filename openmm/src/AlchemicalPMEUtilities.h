#ifndef OPENMM_ALCHEMICAL_PME_UTILITIES_H
#define OPENMM_ALCHEMICAL_PME_UTILITIES_H

#include "AlchemicalPMEForce.h"
#include "openmm/Vec3.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <utility>
#include <vector>

namespace OpenMM {
namespace AlchemicalPME {

constexpr double CoulombConstant = 138.935456;
constexpr double Pi = 3.14159265358979323846;
constexpr double SqrtPi = 1.7724538509055160273;

struct ParticleState {
    double charge;
    double dCharge;
    double sigma;
    double dSigma;
    double epsilon;
    double dEpsilon;
};

struct ParticleDeviceParameter {
    double charge0;
    double chargeDelta;
    double sigma0;
    double sigmaDelta;
    double epsilon0;
    double epsilonDelta;
};

struct ExceptionState {
    int particle1;
    int particle2;
    double chargeProd;
    double dChargeProd;
    double sigma;
    double dSigma;
    double epsilon;
    double dEpsilon;
};

struct ExceptionDeviceParameter {
    int particle1;
    int particle2;
    double chargeProd0;
    double chargeProdDelta;
    double sigma0;
    double sigmaDelta;
    double epsilon0;
    double epsilonDelta;
};

inline double interpolate(double value0, double value1, double lambda)
{
    return value0 + lambda*(value1 - value0);
}

inline double derivative(double value0, double value1)
{
    return value1 - value0;
}

inline double clamp01(double value)
{
    return std::max(0.0, std::min(1.0, value));
}

inline void interpolateElectrostaticEndpoint(double value0,
                                             double value1,
                                             double lambda,
                                             double electrostaticLambdaStart,
                                             double& value,
                                             double& dValue)
{
    if (electrostaticLambdaStart <= 0.0) {
        value = interpolate(value0, value1, lambda);
        dValue = derivative(value0, value1);
        return;
    }

    const double valueTol = 1.0e-12;
    const bool grow = std::abs(value0) <= valueTol && std::abs(value1) > valueTol;
    const bool vanish = std::abs(value1) <= valueTol && std::abs(value0) > valueTol;
    if (grow) {
        const double denom = 1.0 - electrostaticLambdaStart;
        if (lambda <= electrostaticLambdaStart) {
            value = 0.0;
            dValue = 0.0;
        }
        else {
            const double schedule = clamp01((lambda - electrostaticLambdaStart)/denom);
            value = value1*schedule;
            dValue = schedule < 1.0 ? value1/denom : 0.0;
        }
        return;
    }
    if (vanish) {
        const double end = 1.0 - electrostaticLambdaStart;
        if (lambda >= end) {
            value = 0.0;
            dValue = 0.0;
        }
        else {
            const double schedule = clamp01(1.0 - lambda/end);
            value = value0*schedule;
            dValue = schedule > 0.0 ? -value0/end : 0.0;
        }
        return;
    }

    value = interpolate(value0, value1, lambda);
    dValue = derivative(value0, value1);
}

inline double dot(const Vec3& a, const Vec3& b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline Vec3 cross(const Vec3& a, const Vec3& b)
{
    return Vec3(a[1]*b[2] - a[2]*b[1],
                a[2]*b[0] - a[0]*b[2],
                a[0]*b[1] - a[1]*b[0]);
}

inline std::pair<int, int> orderedPair(int particle1, int particle2)
{
    return particle1 < particle2 ? std::make_pair(particle1, particle2) : std::make_pair(particle2, particle1);
}

inline void buildParticleStates(const AlchemicalPMEForce& force,
                                double lambda,
                                std::vector<ParticleState>& particleState)
{
    particleState.resize(force.getNumParticles());
    const double electrostaticLambdaStart = force.getElectrostaticLambdaStart();
    for (int i = 0; i < force.getNumParticles(); ++i) {
        double charge0, charge1, sigma0, sigma1, epsilon0, epsilon1;
        force.getParticleParameters(i, charge0, charge1, sigma0, sigma1, epsilon0, epsilon1);
        double charge, dCharge;
        interpolateElectrostaticEndpoint(charge0, charge1, lambda, electrostaticLambdaStart,
                                         charge, dCharge);
        particleState[i] = {
            charge, dCharge,
            interpolate(sigma0, sigma1, lambda), derivative(sigma0, sigma1),
            interpolate(epsilon0, epsilon1, lambda), derivative(epsilon0, epsilon1)
        };
    }
}

inline void buildExceptionStates(const AlchemicalPMEForce& force,
                                 double lambda,
                                 std::vector<ExceptionState>& exceptionState)
{
    exceptionState.resize(force.getNumExceptions());
    const double electrostaticLambdaStart = force.getElectrostaticLambdaStart();
    for (int i = 0; i < force.getNumExceptions(); ++i) {
        int particle1, particle2;
        double chargeProd0, chargeProd1, sigma0, sigma1, epsilon0, epsilon1;
        force.getExceptionParameters(i, particle1, particle2, chargeProd0, chargeProd1,
                                     sigma0, sigma1, epsilon0, epsilon1);
        double chargeProd, dChargeProd;
        interpolateElectrostaticEndpoint(chargeProd0, chargeProd1, lambda, electrostaticLambdaStart,
                                         chargeProd, dChargeProd);
        exceptionState[i] = {
            particle1, particle2,
            chargeProd, dChargeProd,
            interpolate(sigma0, sigma1, lambda), derivative(sigma0, sigma1),
            interpolate(epsilon0, epsilon1, lambda), derivative(epsilon0, epsilon1)
        };
    }
}

inline bool buildExceptionPairs(const AlchemicalPMEForce& force,
                                std::set<std::pair<int, int>>& exceptionPairs)
{
    exceptionPairs.clear();
    for (int i = 0; i < force.getNumExceptions(); ++i) {
        int particle1, particle2;
        double chargeProd0, chargeProd1, sigma0, sigma1, epsilon0, epsilon1;
        force.getExceptionParameters(i, particle1, particle2, chargeProd0, chargeProd1,
                                     sigma0, sigma1, epsilon0, epsilon1);
        if (!exceptionPairs.insert(orderedPair(particle1, particle2)).second) {
            return false;
        }
    }
    return true;
}

inline void buildParticleDeviceParameters(const AlchemicalPMEForce& force,
                                          std::vector<ParticleDeviceParameter>& parameters)
{
    parameters.resize(force.getNumParticles());
    for (int i = 0; i < force.getNumParticles(); ++i) {
        double charge0, charge1, sigma0, sigma1, epsilon0, epsilon1;
        force.getParticleParameters(i, charge0, charge1, sigma0, sigma1, epsilon0, epsilon1);
        parameters[i] = {
            charge0, charge1 - charge0,
            sigma0, sigma1 - sigma0,
            epsilon0, epsilon1 - epsilon0
        };
    }
}

inline void buildExceptionDeviceParameters(const AlchemicalPMEForce& force,
                                           std::vector<ExceptionDeviceParameter>& parameters)
{
    parameters.resize(force.getNumExceptions());
    for (int i = 0; i < force.getNumExceptions(); ++i) {
        int particle1, particle2;
        double chargeProd0, chargeProd1, sigma0, sigma1, epsilon0, epsilon1;
        force.getExceptionParameters(i, particle1, particle2, chargeProd0, chargeProd1,
                                     sigma0, sigma1, epsilon0, epsilon1);
        parameters[i] = {
            particle1, particle2,
            chargeProd0, chargeProd1 - chargeProd0,
            sigma0, sigma1 - sigma0,
            epsilon0, epsilon1 - epsilon0
        };
    }
}

} // namespace AlchemicalPME
} // namespace OpenMM

#endif
