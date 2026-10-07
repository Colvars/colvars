#include "CudaGlobalMasterClient.h"
#include "colvar_gpu_support.h"
#include "colvarproxy_cudaglobalmaster.h"
#include "colvarproxy_namd.h"
#include "SimParameters.h"
#include "common.h"
#include "ScriptTcl.h"
#include "colvarscript.h"

#include <cstring>

#if defined (NAMD_CUDA) || defined (NAMD_HIP)

#if defined (__linux__) || defined (__APPLE__)
extern "C" {
  CudaGlobalMasterColvars* allocator() {
    return new CudaGlobalMasterColvars();
  }
  void deleter(CudaGlobalMasterColvars* ptr) {
    if (ptr != nullptr) {
      delete ptr;
    }
  }
}
#endif

#ifdef WIN32
extern "C" {
  __declspec (dllexport) CudaGlobalMasterColvars* allocator() {
    return new CudaGlobalMasterColvars();
  }
  __declspec (dllexport) void deleter(CudaGlobalMasterColvars* ptr) {
    if (ptr != nullptr) {
      delete ptr;
    }
  }
}
#endif

CudaGlobalMasterColvars::CudaGlobalMasterColvars():
  CudaGlobalMasterClient()
{
  if (CudaGlobalMasterClient::getSimParameters()->colvarsOn) {
    NAMD_die("This plugin is incompatible with the Colvars bundled with NAMD.");
  }
  mImpl = std::make_unique<colvarproxy_namd>();
}

CudaGlobalMasterColvars::~CudaGlobalMasterColvars() {}

void CudaGlobalMasterColvars::initialize(
  const std::vector<std::string>& arguments,
  int deviceID, cudaStream_t stream) {
  CudaGlobalMasterClient::initialize(arguments, deviceID, stream);
  mDeviceID = deviceID;
  const int error_code = mImpl->initialize_from_cudagm(
    this, arguments, deviceID, stream,
    CudaGlobalMasterClient::getSimParameters(),
    CudaGlobalMasterClient::getMolecule(),
    CudaGlobalMasterClient::getScript());
  if (error_code != COLVARS_OK) {
    NAMD_die("Error initializing the CUDA GlobalMaster Colvars proxy.");
  }
}

bool CudaGlobalMasterColvars::requestedAtomsChanged()  {
  return mImpl->atomsChanged();
}

const std::vector<AtomID>& CudaGlobalMasterColvars::getRequestedTotalForcesAtoms() const {
  if (mImpl->total_forces_enabled()) {
    return *(mImpl->get_atom_ids());
  } else {
    return mEmpty;
  }
}

bool CudaGlobalMasterColvars::requestedTotalForcesAtomsChanged() {
  // The server must also be notified when total-force requests are removed.
  return mImpl->atomsChanged();
}

bool CudaGlobalMasterColvars::requestUpdateMasses() {
  return mImpl->atomsChanged();
}

bool CudaGlobalMasterColvars::requestUpdateCharges() {
  return mImpl->atomsChanged();
}

#if CUDAGM_VERSION >= 3
void CudaGlobalMasterColvars::setStep(int64_t step, int startup, int doMigration) {
  CudaGlobalMasterClient::setStep(step, startup, doMigration);
#else
void CudaGlobalMasterColvars::setStep(int64_t step) {
  CudaGlobalMasterClient::setStep(step);
#endif
  if (mImpl->atomsChanged()) {
    mImpl->reallocate();
  }
}

void CudaGlobalMasterColvars::calculate() {
  mImpl->calculate();
}

cudaStream_t CudaGlobalMasterColvars::getStream() {
  return mImpl->getStream();
}

bool CudaGlobalMasterColvars::requestUpdateAtomTotalForces() {
  return mImpl->total_forces_enabled();
}

double CudaGlobalMasterColvars::getEnergy() const {
  return mImpl->getEnergy();
}

double* CudaGlobalMasterColvars::getAppliedForces() const {
  return mImpl->getAppliedForces();
}

double* CudaGlobalMasterColvars::getPositions() {
  return mImpl->getPositions();
}

float* CudaGlobalMasterColvars::getMasses() {
  return mImpl->getMasses();
}

float* CudaGlobalMasterColvars::getCharges() {
  return mImpl->getCharges();
}

double* CudaGlobalMasterColvars::getTotalForces() {
  return mImpl->getTotalForces();
}

double* CudaGlobalMasterColvars::getLattice() {
  return mImpl->getLattice();
}

const std::vector<AtomID>& CudaGlobalMasterColvars::getRequestedAtoms() const {
  return *(mImpl->get_atom_ids());
}

void CudaGlobalMasterColvars::onBuffersUpdated() {
  mImpl->onBuffersUpdated();
}

int CudaGlobalMasterColvars::updateFromTCLCommand(const std::vector<std::string>& arguments) {
  mTCLResult.clear();
  if (mDeviceID < 0 || arguments.size() < 2) {
    return TCL_ERROR;
  }

  // The proxy's guard is file-local. Keep the client's device selected for the
  // entire script command, which can create or destroy GPU-backed components.
  class device_guard {
  public:
    explicit device_guard(int device) {
      if (colvars_gpu::gpuAssert(cudaGetDevice(&savedDevice), __FILE__, __LINE__) != COLVARS_OK ||
          colvars_gpu::gpuAssert(cudaSetDevice(device), __FILE__, __LINE__) != COLVARS_OK) {
        NAMD_die("Cannot select the Colvars CUDA GlobalMaster scripting device.");
      }
    }
    ~device_guard() {
      if (colvars_gpu::gpuAssert(cudaSetDevice(savedDevice), __FILE__, __LINE__) != COLVARS_OK) {
        NAMD_die("Cannot restore the CUDA device after a Colvars script command.");
      }
    }
    device_guard(device_guard const &) = delete;
    device_guard &operator=(device_guard const &) = delete;
  private:
    int savedDevice = -1;
  } guard(mDeviceID);

  // Prepare the arguments
  const int objc = arguments.size();
  unsigned char** objv = new unsigned char*[objc];
  for (int i = 0; i < objc; ++i) {
    const int len = std::strlen(arguments[i].c_str());
    objv[i] = new unsigned char[len + 1];
    std::strncpy(reinterpret_cast<char*>(objv[i]),
                 arguments[i].c_str(), len + 1);
    objv[i][len] = '\0';
  }
  // Call Colvars scripting interface
  int error_code = TCL_ERROR;
  if (mImpl->script) {
    // objv[0] is the name of this client
    if (COLVARSCRIPT_OK == mImpl->script->run(objc-1, objv+1)) {
      error_code = TCL_OK;
      mTCLResult = mImpl->get_error_msgs() + mImpl->script->str_result();
    }
  }
  // Cleanups
  for (int i = 0; i < objc; ++i) {
    delete[] objv[i];
  }
  delete[] objv;
  return error_code;
}

#endif // defined (NAMD_CUDA) || defined (NAMD_HIP)
