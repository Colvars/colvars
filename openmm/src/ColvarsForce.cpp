#include "ColvarsForce.h"
#include "ColvarsForceImpl.h"

#include "openmm/OpenMMException.h"

#include <cmath>

using namespace OpenMM;

namespace {

void validate_nonempty(const std::string& value, const std::string& name)
{
    if (value.empty()) {
        throw OpenMMException(name + " must not be empty");
    }
}

} // namespace

ColvarsForce::ColvarsForce(const std::string& configFile)
    : ColvarsForce(configFile, "colvars", "lambda")
{
}

ColvarsForce::ColvarsForce(const std::string& configFile,
                           const std::string& outputPrefix,
                           const std::string& lambdaParameterName)
    : configSource_(ConfigSource::File),
      config_(configFile),
      outputPrefix_(outputPrefix),
      restartOutputPrefix_(outputPrefix),
      restartInputPrefix_(),
      lambdaParameterName_(lambdaParameterName),
      lambdaDerivativeParameterName_("dE_dlambda"),
      tabfAlphaParameterName_("tabf_alpha"),
      tabfAlphaDerivativeParameterName_("tabf_dE_dalpha"),
      alchemyEnabled_(true),
      tabfEnabled_(false),
      scaleLambdaDerivativeByTabfAlpha_(false),
      pbcMode_(PbcMode::Auto),
      targetTemperatureK_(300.0),
      randomSeed_(0),
      restartFrequency_(0)
{
    validate_nonempty(config_, "configFile");
    validate_nonempty(outputPrefix_, "outputPrefix");
    validate_nonempty(lambdaParameterName_, "lambdaParameterName");
}

ColvarsForce ColvarsForce::fromConfigString(const std::string& configString,
                                            const std::string& outputPrefix,
                                            const std::string& lambdaParameterName)
{
    ColvarsForce force("unused.colvars", outputPrefix, lambdaParameterName);
    force.setConfigString(configString);
    return force;
}

ColvarsForce::ConfigSource ColvarsForce::getConfigSource() const
{
    return configSource_;
}

const std::string& ColvarsForce::getConfig() const
{
    return config_;
}

void ColvarsForce::setConfigFile(const std::string& configFile)
{
    validate_nonempty(configFile, "configFile");
    configSource_ = ConfigSource::File;
    config_ = configFile;
}

void ColvarsForce::setConfigString(const std::string& configString)
{
    validate_nonempty(configString, "configString");
    configSource_ = ConfigSource::String;
    config_ = configString;
}

const std::string& ColvarsForce::getOutputPrefix() const
{
    return outputPrefix_;
}

void ColvarsForce::setOutputPrefix(const std::string& outputPrefix)
{
    validate_nonempty(outputPrefix, "outputPrefix");
    outputPrefix_ = outputPrefix;
}

const std::string& ColvarsForce::getRestartOutputPrefix() const
{
    return restartOutputPrefix_;
}

void ColvarsForce::setRestartOutputPrefix(const std::string& restartOutputPrefix)
{
    validate_nonempty(restartOutputPrefix, "restartOutputPrefix");
    restartOutputPrefix_ = restartOutputPrefix;
}

const std::string& ColvarsForce::getRestartInputPrefix() const
{
    return restartInputPrefix_;
}

void ColvarsForce::setRestartInputPrefix(const std::string& restartInputPrefix)
{
    restartInputPrefix_ = restartInputPrefix;
}

int ColvarsForce::getRestartFrequency() const
{
    return restartFrequency_;
}

void ColvarsForce::setRestartFrequency(int restartFrequency)
{
    if (restartFrequency < 0) {
        throw OpenMMException("restartFrequency must be non-negative");
    }
    restartFrequency_ = restartFrequency;
}

const std::string& ColvarsForce::getLambdaParameterName() const
{
    return lambdaParameterName_;
}

void ColvarsForce::setLambdaParameterName(const std::string& lambdaParameterName)
{
    validate_nonempty(lambdaParameterName, "lambdaParameterName");
    lambdaParameterName_ = lambdaParameterName;
}

const std::string& ColvarsForce::getLambdaDerivativeParameterName() const
{
    return lambdaDerivativeParameterName_;
}

void ColvarsForce::setLambdaDerivativeParameterName(const std::string& derivativeParameterName)
{
    validate_nonempty(derivativeParameterName, "derivativeParameterName");
    lambdaDerivativeParameterName_ = derivativeParameterName;
}

const std::string& ColvarsForce::getTabfAlphaParameterName() const
{
    return tabfAlphaParameterName_;
}

void ColvarsForce::setTabfAlphaParameterName(const std::string& alphaParameterName)
{
    validate_nonempty(alphaParameterName, "alphaParameterName");
    tabfAlphaParameterName_ = alphaParameterName;
}

const std::string& ColvarsForce::getTabfAlphaDerivativeParameterName() const
{
    return tabfAlphaDerivativeParameterName_;
}

void ColvarsForce::setTabfAlphaDerivativeParameterName(const std::string& derivativeParameterName)
{
    validate_nonempty(derivativeParameterName, "derivativeParameterName");
    tabfAlphaDerivativeParameterName_ = derivativeParameterName;
}

bool ColvarsForce::getTabfEnabled() const
{
    return tabfEnabled_;
}

void ColvarsForce::setTabfEnabled(bool enabled)
{
    tabfEnabled_ = enabled;
}

bool ColvarsForce::getScaleLambdaDerivativeByTabfAlpha() const
{
    return scaleLambdaDerivativeByTabfAlpha_;
}

void ColvarsForce::setScaleLambdaDerivativeByTabfAlpha(bool enabled)
{
    scaleLambdaDerivativeByTabfAlpha_ = enabled;
}

bool ColvarsForce::getAlchemyEnabled() const
{
    return alchemyEnabled_;
}

void ColvarsForce::setAlchemyEnabled(bool enabled)
{
    alchemyEnabled_ = enabled;
}

double ColvarsForce::getTargetTemperature() const
{
    return targetTemperatureK_;
}

void ColvarsForce::setTargetTemperature(double temperatureK)
{
    if (!std::isfinite(temperatureK) || temperatureK < 0.0) {
        throw OpenMMException("temperatureK must be finite and non-negative");
    }
    targetTemperatureK_ = temperatureK;
}

int ColvarsForce::getRandomSeed() const
{
    return randomSeed_;
}

void ColvarsForce::setRandomSeed(int seed)
{
    if (seed < 0) {
        throw OpenMMException("random seed must be non-negative");
    }
    randomSeed_ = seed;
}

const std::string& ColvarsForce::getUnitSystem() const
{
    static const std::string unitSystem{"gromacs"};
    return unitSystem;
}

int ColvarsForce::getNumParticleCharges() const
{
    return static_cast<int>(particleCharges_.size());
}

double ColvarsForce::getParticleCharge(int particleIndex) const
{
    if (particleIndex < 0 || particleIndex >= static_cast<int>(particleCharges_.size())) {
        return 0.0;
    }
    return particleCharges_[particleIndex];
}

void ColvarsForce::setParticleCharge(int particleIndex, double chargeElementary)
{
    if (particleIndex < 0) {
        throw OpenMMException("particleIndex must be non-negative");
    }
    if (!std::isfinite(chargeElementary)) {
        throw OpenMMException("chargeElementary must be finite");
    }
    if (particleIndex >= static_cast<int>(particleCharges_.size())) {
        particleCharges_.resize(static_cast<size_t>(particleIndex) + 1, 0.0);
    }
    particleCharges_[particleIndex] = chargeElementary;
}

void ColvarsForce::setParticleCharges(const std::vector<double>& chargesElementary)
{
    for (double charge : chargesElementary) {
        if (!std::isfinite(charge)) {
            throw OpenMMException("particle charges must be finite");
        }
    }
    particleCharges_ = chargesElementary;
}

bool ColvarsForce::usesPeriodicBoundaryConditions() const
{
    return pbcMode_ == PbcMode::Enabled;
}

void ColvarsForce::setUsesPeriodicBoundaryConditions(bool usesPbc)
{
    pbcMode_ = usesPbc ? PbcMode::Enabled : PbcMode::Disabled;
}

ColvarsForce::PbcMode ColvarsForce::getPbcMode() const
{
    return pbcMode_;
}

void ColvarsForce::setPbcMode(PbcMode mode)
{
    switch (mode) {
        case PbcMode::Auto:
        case PbcMode::Enabled:
        case PbcMode::Disabled:
            pbcMode_ = mode;
            return;
    }
    throw OpenMMException("Invalid ColvarsForce PBC mode");
}

ForceImpl* ColvarsForce::createImpl() const
{
    return new ColvarsForceImpl(*this);
}
