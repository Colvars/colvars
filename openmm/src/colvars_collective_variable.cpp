#include "ColvarsCollectiveVariable.h"
#include "ColvarsCollectiveVariableImpl.h"
#include "colvarproxy_openmm.h"

#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/serialization/SerializationProxy.h"

#include <algorithm>
#include <regex>
#include <typeinfo>

using namespace OpenMM;

namespace {

void validate_nonempty(const std::string& value, const std::string& name)
{
    if (value.empty()) {
        throw OpenMMException("ColvarsCollectiveVariable: " + name + " must not be empty");
    }
}

// Silences most of the Colvars log unless the configuration sets logLevel, which
// Colvars does not allow twice
std::string with_default_log_level(const std::string& config)
{
    static const std::regex logLevel(R"((^|\n)[ \t]*loglevel\b)", std::regex::icase);
    return std::regex_search(config, logLevel) ? config : "logLevel 0\n" + config;
}

ColvarsForce make_settings(const ColvarsCollectiveVariable& owner)
{
    // The output prefix is required by ColvarsForce but unused: the proxy writes no
    // files when it reports a colvar
    ColvarsForce settings = ColvarsForce::fromConfigString(
        with_default_log_level(owner.getConfigString()), "colvars");
    settings.setAlchemyEnabled(false);
    settings.setUsesPeriodicBoundaryConditions(owner.usesPeriodicBoundaryConditions());
    return settings;
}

class ColvarsCollectiveVariableProxy final : public SerializationProxy {
public:
    ColvarsCollectiveVariableProxy() : SerializationProxy("ColvarsCollectiveVariable") {}

    void serialize(const void* object, SerializationNode& node) const override
    {
        const auto& force = *reinterpret_cast<const ColvarsCollectiveVariable*>(object);
        node.setIntProperty("version", 1);
        node.setIntProperty("forceGroup", force.getForceGroup());
        node.setStringProperty("name", force.getName());
        node.setStringProperty("configString", force.getConfigString());
        node.setStringProperty("colvarName", force.getColvarName());
        node.setBoolProperty("usesPeriodic", force.usesPeriodicBoundaryConditions());
    }

    void* deserialize(const SerializationNode& node) const override
    {
        if (node.getIntProperty("version") != 1) {
            throw OpenMMException("Unsupported version number");
        }
        auto* force = new ColvarsCollectiveVariable(node.getStringProperty("configString"),
                                                    node.getStringProperty("colvarName"));
        force->setForceGroup(node.getIntProperty("forceGroup", 0));
        force->setName(node.getStringProperty("name", force->getName()));
        force->setUsesPeriodicBoundaryConditions(node.getBoolProperty("usesPeriodic", false));
        return force;
    }
};

// Registers the serialization proxy whenever this class is linked in
[[maybe_unused]] const bool proxyRegistered = [] {
    SerializationProxy::registerProxy(typeid(ColvarsCollectiveVariable),
                                      new ColvarsCollectiveVariableProxy());
    return true;
}();

} // namespace

ColvarsCollectiveVariable::ColvarsCollectiveVariable(const std::string& configString,
                                                     const std::string& colvarName)
    : configString_(configString), colvarName_(colvarName), usesPbc_(false)
{
    validate_nonempty(configString_, "configString");
    validate_nonempty(colvarName_, "colvarName");
    setName("ColvarsCollectiveVariable");
}

const std::string& ColvarsCollectiveVariable::getConfigString() const
{
    return configString_;
}

const std::string& ColvarsCollectiveVariable::getColvarName() const
{
    return colvarName_;
}

bool ColvarsCollectiveVariable::usesPeriodicBoundaryConditions() const
{
    return usesPbc_;
}

void ColvarsCollectiveVariable::setUsesPeriodicBoundaryConditions(bool usesPbc)
{
    usesPbc_ = usesPbc;
}

ForceImpl* ColvarsCollectiveVariable::createImpl() const
{
    return new ColvarsCollectiveVariableImpl(*this);
}

ColvarsCollectiveVariableImpl::ColvarsCollectiveVariableImpl(
    const ColvarsCollectiveVariable& owner)
    : CustomCPPForceImpl(owner),
      owner_(owner),
      settings_(make_settings(owner))
{
}

ColvarsCollectiveVariableImpl::~ColvarsCollectiveVariableImpl() = default;

const Force& ColvarsCollectiveVariableImpl::getOwner() const
{
    return owner_;
}

void ColvarsCollectiveVariableImpl::initialize(ContextImpl& context)
{
    if (owner_.getForceGroup() == 31) {
        throw OpenMMException(
            "ColvarsCollectiveVariable group 31 is unsupported because this OpenMM "
            "dependency uses a signed force-group shift. Use group 0-30.");
    }
    CustomCPPForceImpl::initialize(context);
    proxy_ = std::make_unique<colvarproxy_openmm>(
        settings_, context.getSystem(), owner_.getColvarName());
}

double ColvarsCollectiveVariableImpl::computeForce(ContextImpl& context,
                                                   const std::vector<Vec3>& positions,
                                                   std::vector<Vec3>& forces)
{
    std::fill(forces.begin(), forces.end(), Vec3());
    const double energy = proxy_->calculate(context, positions);
    proxy_->addForcesToOpenMM(forces);
    return energy;
}

std::map<std::string, double> ColvarsCollectiveVariableImpl::getDefaultParameters()
{
    return {};
}
