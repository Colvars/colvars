%module openmm_colvars

%pythonbegin %{
import os as _openmm_colvars_os
import sys as _openmm_colvars_sys

if _openmm_colvars_sys.platform == "win32" and hasattr(_openmm_colvars_os, "add_dll_directory"):
    for _openmm_colvars_dll_relpath in (
            ("Library", "lib"),
            ("Library", "bin"),
            ("Library", "lib", "plugins")):
        _openmm_colvars_dll_dir = _openmm_colvars_os.path.join(
            _openmm_colvars_sys.prefix, *_openmm_colvars_dll_relpath)
        if _openmm_colvars_os.path.isdir(_openmm_colvars_dll_dir):
            _openmm_colvars_os.add_dll_directory(_openmm_colvars_dll_dir)

try:
    import openmm as _openmm_colvars_openmm
except Exception:
    _openmm_colvars_openmm = None
%}

%include "std_string.i"
%include "std_vector.i"
%include "typemaps.i"

%exception {
    try {
        $action
    } catch (const std::exception& exception) {
        OpenMMColvars_SetPythonException(exception.what());
        SWIG_fail;
    } catch (...) {
        PyErr_SetString(PyExc_RuntimeError, "Unknown C++ exception in openmm_colvars");
        SWIG_fail;
    }
}

%apply int &OUTPUT { int& particle1, int& particle2 };
%apply double &OUTPUT {
    double& charge0, double& charge1,
    double& sigma0, double& sigma1,
    double& epsilon0, double& epsilon1,
    double& chargeProd0, double& chargeProd1
};

%{
#define SWIG_FILE_WITH_INIT
#include <exception>
#include "AlchemicalPMEForce.h"
#include "ColvarsCollectiveVariable.h"
#include "openmm/serialization/XmlSerializer.h"
#include <memory>
#include <sstream>
#include "ColvarsForce.h"
#include "ColvarsTABFIntegrator.h"

static void OpenMMColvars_SetPythonException(const char* message) {
    PyObject* module = PyImport_ImportModule("openmm");
    PyObject* exceptionType = module == nullptr ? nullptr :
        PyObject_GetAttrString(module, "OpenMMException");
    if (exceptionType != nullptr && PyExceptionClass_Check(exceptionType)) {
        PyErr_SetString(exceptionType, message);
    } else {
        PyErr_Clear();
        PyErr_SetString(PyExc_RuntimeError, message);
    }
    Py_XDECREF(exceptionType);
    Py_XDECREF(module);
}
%}

%import(module="openmm.openmm") "openmm_colvars_openmm_import.i"

namespace std {
    %template(vector_double) vector<double>;
}

namespace ColvarsPlugin {

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
    void setPMEParameters(double alpha, int nx, int ny, int nz);
    int getPMEInterpolationOrder() const;
    void setPMEInterpolationOrder(int order);
    bool usesPeriodicBoundaryConditions() const;
};

class ColvarsForce : public OpenMM::Force {
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

    bool usesPeriodicBoundaryConditions() const;
    void setUsesPeriodicBoundaryConditions(bool usesPbc);
    PbcMode getPbcMode() const;
    void setPbcMode(PbcMode mode);
};

class ColvarsCollectiveVariable : public Force {
public:
    ColvarsCollectiveVariable(const std::string& configString, const std::string& colvarName);
    ColvarsCollectiveVariable(const ColvarsCollectiveVariable& other);

    const std::string& getConfigString() const;
    const std::string& getColvarName() const;

    bool usesPeriodicBoundaryConditions() const;
    void setUsesPeriodicBoundaryConditions(bool usesPbc);
};

%newobject OpenMM::ColvarsCollectiveVariable::_deserialize;
%extend ColvarsCollectiveVariable {
    // XmlSerializer.deserialize returns a plain Force wrapper for plugin types, so
    // unpickling needs a typed object
    static ColvarsCollectiveVariable* _deserialize(const std::string& xml) {
        std::istringstream stream(xml);
        std::unique_ptr<OpenMM::Force> force(OpenMM::XmlSerializer::deserialize<OpenMM::Force>(stream));
        auto* cv = dynamic_cast<OpenMM::ColvarsCollectiveVariable*>(force.get());
        if (cv == nullptr) {
            throw OpenMM::OpenMMException("The XML does not describe a ColvarsCollectiveVariable");
        }
        force.release();
        return cv;
    }

    %pythoncode %{
    def __copy__(self):
        return ColvarsCollectiveVariable(self)

    def __deepcopy__(self, memo):
        return ColvarsCollectiveVariable(self)

    def __getstate__(self):
        return _openmm_colvars_openmm.XmlSerializer.serialize(self)

    def __setstate__(self, state):
        self.this = ColvarsCollectiveVariable._deserialize(state).this
    %}
}

class ColvarsTABFIntegrator : public OpenMM::CustomIntegrator {
public:
    ColvarsTABFIntegrator(double stepSize,
                          double temperature,
                          double frictionCoeff,
                          double initialAlpha = 1.0,
                          int physicalForceGroup = 0,
                          int colvarsForceGroup = 1,
                          const std::string& alphaVariableName = "tabf_alpha",
                          const std::string& alphaDerivativeVariableName = "tabf_dE_dalpha");

    int getPhysicalForceGroup() const;
    int getColvarsForceGroup() const;
    const std::string& getAlphaVariableName() const;
    const std::string& getAlphaDerivativeVariableName() const;
};

} // namespace ColvarsPlugin
