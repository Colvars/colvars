%module openmm

%{
#include "openmm/Force.h"
#include "openmm/Integrator.h"
#include "openmm/CustomIntegrator.h"
%}

namespace OpenMM {

class Force {
public:
    virtual ~Force();
    int getForceGroup() const;
    void setForceGroup(int group);
    virtual bool usesPeriodicBoundaryConditions() const;
};

class Integrator {
public:
    virtual ~Integrator();
};

class CustomIntegrator : public Integrator {
public:
    CustomIntegrator(double stepSize);
    int addGlobalVariable(const std::string& name, double initialValue);
    double getGlobalVariableByName(const std::string& name) const;
    void setGlobalVariableByName(const std::string& name, double value);
};

} // namespace OpenMM
