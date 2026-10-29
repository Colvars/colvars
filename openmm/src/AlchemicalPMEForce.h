#ifndef OPENMM_ALCHEMICAL_PME_FORCE_H
#define OPENMM_ALCHEMICAL_PME_FORCE_H

#include "openmm/Force.h"

#include <string>
#include <vector>

namespace ColvarsPlugin {

class AlchemicalPMEForceImpl;

/**
 * PME nonbonded force with particle parameters linearly interpolated by a
 * global alchemical parameter.
 *
 * This first implementation is intentionally conservative: it exposes a small
 * NAMD/OpenMM-FEP style parameter model and computes an analytic dU/dlambda
 * companion context parameter for Colvars alchLambda/eABF.
 */
class AlchemicalPMEForce : public OpenMM::Force {
public:
    enum NonbondedMethod {
        NoCutoff = 0,
        CutoffPeriodic = 1,
        PME = 2
    };

    AlchemicalPMEForce();

    int addParticle(double charge0, double charge1,
                    double sigma0, double sigma1,
                    double epsilon0, double epsilon1);

    int getNumParticles() const;
    void getParticleParameters(int index,
                               double& charge0, double& charge1,
                               double& sigma0, double& sigma1,
                               double& epsilon0, double& epsilon1) const;
    void setParticleParameters(int index,
                               double charge0, double charge1,
                               double sigma0, double sigma1,
                               double epsilon0, double epsilon1);

    int addException(int particle1, int particle2,
                     double chargeProd0, double chargeProd1,
                     double sigma0, double sigma1,
                     double epsilon0, double epsilon1);
    int getNumExceptions() const;
    void getExceptionParameters(int index,
                                int& particle1, int& particle2,
                                double& chargeProd0, double& chargeProd1,
                                double& sigma0, double& sigma1,
                                double& epsilon0, double& epsilon1) const;
    void setExceptionParameters(int index,
                                int particle1, int particle2,
                                double chargeProd0, double chargeProd1,
                                double sigma0, double sigma1,
                                double epsilon0, double epsilon1);

    const std::string& getLambdaParameterName() const;
    void setLambdaParameterName(const std::string& name);

    const std::string& getLambdaDerivativeParameterName() const;
    void setLambdaDerivativeParameterName(const std::string& name);

    double getElectrostaticLambdaStart() const;
    void setElectrostaticLambdaStart(double lambdaStart);

    NonbondedMethod getNonbondedMethod() const;
    void setNonbondedMethod(NonbondedMethod method);

    double getCutoffDistance() const;
    void setCutoffDistance(double distance);
    double getReactionFieldDielectric() const;
    void setReactionFieldDielectric(double dielectric);

    void getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const;
    void setPMEParameters(double alpha, int nx, int ny, int nz);
    int getPMEInterpolationOrder() const;
    void setPMEInterpolationOrder(int order);

    bool usesPeriodicBoundaryConditions() const override;

protected:
    OpenMM::ForceImpl* createImpl() const override;

private:
    struct ParticleInfo {
        double charge0;
        double charge1;
        double sigma0;
        double sigma1;
        double epsilon0;
        double epsilon1;
    };

    struct ExceptionInfo {
        int particle1;
        int particle2;
        double chargeProd0;
        double chargeProd1;
        double sigma0;
        double sigma1;
        double epsilon0;
        double epsilon1;
    };

    std::vector<ParticleInfo> particles_;
    std::vector<ExceptionInfo> exceptions_;
    std::string lambdaParameterName_;
    std::string lambdaDerivativeParameterName_;
    double electrostaticLambdaStart_;
    NonbondedMethod nonbondedMethod_;
    double cutoffDistance_;
    double reactionFieldDielectric_;
    double pmeAlpha_;
    int pmeGridX_;
    int pmeGridY_;
    int pmeGridZ_;
    int pmeInterpolationOrder_;
};

} // namespace ColvarsPlugin

#endif
