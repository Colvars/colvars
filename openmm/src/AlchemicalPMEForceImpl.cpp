#include "AlchemicalPMEForceImpl.h"

#include "AlchemicalPMEKernels.h"
#include "openmm/internal/ContextImpl.h"

#include <cstdint>

using namespace OpenMM;
using namespace ColvarsPlugin;

AlchemicalPMEForceImpl::AlchemicalPMEForceImpl(const AlchemicalPMEForce& owner)
    : owner_(owner)
{
    forceGroup = owner.getForceGroup();
}

const AlchemicalPMEForce& AlchemicalPMEForceImpl::getOwner() const
{
    return owner_;
}

void AlchemicalPMEForceImpl::initialize(ContextImpl& context)
{
    registerAlchemicalPMEKernelFactories();
    kernel_ = context.getPlatform().createKernel(CalcAlchemicalPMEForceKernel::Name(), context);
    kernel_.getAs<CalcAlchemicalPMEForceKernel>().initialize(context, owner_);
}

void AlchemicalPMEForceImpl::updateContextState(ContextImpl& context, bool& forcesInvalid)
{
}

double AlchemicalPMEForceImpl::calcForcesAndEnergy(ContextImpl& context,
                                                   bool includeForces,
                                                   bool includeEnergy,
                                                   int groups)
{
    const std::uint32_t forceGroupFlag =
        std::uint32_t{1} << static_cast<unsigned>(forceGroup);
    if ((static_cast<std::uint32_t>(groups) & forceGroupFlag) == 0) {
        return 0.0;
    }
    return kernel_.getAs<CalcAlchemicalPMEForceKernel>().execute(context, includeForces, includeEnergy);
}

std::map<std::string, double> AlchemicalPMEForceImpl::getDefaultParameters()
{
    return {
        {owner_.getLambdaParameterName(), 0.0}
    };
}

std::vector<std::string> AlchemicalPMEForceImpl::getKernelNames()
{
    registerAlchemicalPMEKernelFactories();
    return {CalcAlchemicalPMEForceKernel::Name()};
}
