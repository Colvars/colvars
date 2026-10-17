#ifndef OPENMM_COLVARS_FORCE_H
#define OPENMM_COLVARS_FORCE_H

#include "openmm/Force.h"

#include <string>
#include <vector>

namespace OpenMM {

class ColvarsForceImpl;

/// OpenMM Force wrapper for the Colvars collective variable module.
///
/// The Colvars configuration uses the native OpenMM/GROMACS unit system by
/// default: nanometers for length and kJ/mol for energy. OpenMM positions,
/// energies, and forces are passed through without a custom conversion layer.
class ColvarsForce : public Force {
public:
    enum class ConfigSource {
        File,
        String
    };

    enum class PbcMode {
        Auto,
        Enabled,
        Disabled
    };

    explicit ColvarsForce(const std::string& configFile);

    ColvarsForce(const std::string& configFile,
                 const std::string& outputPrefix,
                 const std::string& lambdaParameterName = "lambda");

    static ColvarsForce fromConfigString(const std::string& configString,
                                         const std::string& outputPrefix,
                                         const std::string& lambdaParameterName = "lambda");

    ConfigSource getConfigSource() const;
    const std::string& getConfig() const;

    void setConfigFile(const std::string& configFile);
    void setConfigString(const std::string& configString);

    const std::string& getOutputPrefix() const;
    void setOutputPrefix(const std::string& outputPrefix);

    const std::string& getRestartOutputPrefix() const;
    void setRestartOutputPrefix(const std::string& restartOutputPrefix);

    const std::string& getRestartInputPrefix() const;
    void setRestartInputPrefix(const std::string& restartInputPrefix);

    int getRestartFrequency() const;
    void setRestartFrequency(int restartFrequency);

    const std::string& getLambdaParameterName() const;
    void setLambdaParameterName(const std::string& lambdaParameterName);

    const std::string& getLambdaDerivativeParameterName() const;
    void setLambdaDerivativeParameterName(const std::string& derivativeParameterName);

    const std::string& getTabfAlphaParameterName() const;
    void setTabfAlphaParameterName(const std::string& alphaParameterName);

    const std::string& getTabfAlphaDerivativeParameterName() const;
    void setTabfAlphaDerivativeParameterName(const std::string& derivativeParameterName);

    bool getTabfEnabled() const;
    void setTabfEnabled(bool enabled);

    bool getScaleLambdaDerivativeByTabfAlpha() const;
    void setScaleLambdaDerivativeByTabfAlpha(bool enabled);

    bool getAlchemyEnabled() const;
    void setAlchemyEnabled(bool enabled);

    double getTargetTemperature() const;
    void setTargetTemperature(double temperatureK);

    int getRandomSeed() const;
    void setRandomSeed(int seed);

    const std::string& getUnitSystem() const;

    int getNumParticleCharges() const;
    double getParticleCharge(int particleIndex) const;
    void setParticleCharge(int particleIndex, double chargeElementary);
    void setParticleCharges(const std::vector<double>& chargesElementary);

    bool usesPeriodicBoundaryConditions() const override;
    void setUsesPeriodicBoundaryConditions(bool usesPbc);
    PbcMode getPbcMode() const;
    void setPbcMode(PbcMode mode);

protected:
    ForceImpl* createImpl() const override;

private:
    ConfigSource configSource_;
    std::string config_;
    std::string outputPrefix_;
    std::string restartOutputPrefix_;
    std::string restartInputPrefix_;
    std::string lambdaParameterName_;
    std::string lambdaDerivativeParameterName_;
    std::string tabfAlphaParameterName_;
    std::string tabfAlphaDerivativeParameterName_;
    bool alchemyEnabled_;
    bool tabfEnabled_;
    bool scaleLambdaDerivativeByTabfAlpha_;
    PbcMode pbcMode_;
    double targetTemperatureK_;
    int randomSeed_;
    int restartFrequency_;
    std::vector<double> particleCharges_;
};

} // namespace OpenMM

#endif
