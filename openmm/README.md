# OpenMM Colvars Interface

OpenMM Colvars Interface adds the [Colvars](https://github.com/Colvars/colvars)
collective-variable module to OpenMM through a C++ `Force` implementation and
Python bindings. It also includes `openmm_fep`, a Python API for NAMD-compatible
alchemical free-energy workflows.

The project supports OpenMM's CPU and CUDA platforms. CUDA builds can enable a
deferred device path for ABF-family total-force sampling.

## Contents

- [Features](#features)
- [Platform support](#platform-support)
- [Requirements](#requirements)
- [Installation in an existing conda environment](#installation-in-an-existing-conda-environment)
- [Quick verification](#quick-verification)
- [Usage](#usage)
- [OpenMM FEP helpers](#openmm-fep-helpers)
- [Units](#units)
- [Limitations](#limitations)

## Features

- Geometry-based Colvars operating on OpenMM particle positions.
- Colvars bias energies and forces applied back to OpenMM.
- ABF-family total-force collection.
- TABF `tabfAlpha` communication when supported by the selected Colvars source.
- `ColvarsTABFIntegrator` for OpenMM-side `U/alpha` dynamics.
- Orthorhombic and triclinic periodic boxes.
- Alchemical lambda through OpenMM context parameters.
- Colvars `customFunction` expressions through OpenMM's Lepton library.
- NAMD-style dual-topology FEP helpers and `.fepout` output.
- C++ and Python interfaces.

## Platform support

| OpenMM platform | Geometry biases | ABF/eABF total forces | Build setting |
| --- | --- | --- | --- |
| CPU | Supported | Supported | `OPENMM_COLVARS_ENABLE_CUDA_KERNELS=OFF` |
| Reference | Supported | Supported | `OPENMM_COLVARS_ENABLE_CUDA_KERNELS=OFF` |
| CUDA | Supported | Supported with the deferred CUDA path | `OPENMM_COLVARS_ENABLE_CUDA_KERNELS=ON` |
| OpenCL / HIP | Supported for ordinary geometry biases | Not supported for same-step total-force sampling | CUDA path unavailable |

For CPU and Reference, add `ColvarsForce` after all physical forces that should
contribute to ABF, ideally as the final force in the `System`. CUDA ABF/eABF
requires both OpenMM's CUDA internal headers and the `OpenMMCUDA` library from
the same OpenMM installation.

Ordinary geometry biases do not require the deferred CUDA path. On platforms
without that path, use force groups 0-30 for `ColvarsForce`; force group 31 is
not supported by the upstream `CustomCPPForceImpl`.

## Requirements

- Linux or WSL with a Bash shell for the provided installation script.
- An activated conda or micromamba environment containing OpenMM.
- A C++17 compiler, CMake, SWIG, and Make.
- OpenMM and Colvars source trees.
- For CUDA builds: an OpenMM installation with its CUDA plugin, NVIDIA CUDA
  driver/development headers (`cuda.h` and `builtin_types.h`), and a compatible
  NVIDIA driver.

The OpenMM source version should match the OpenMM package installed in the
active environment. Use a recent stable Colvars source tree; Colvars treats its
`master` branch as the stable, continuously updated branch.

## Installation in an existing conda environment

This is the recommended installation method. The script links the extension
against the OpenMM already installed in the active environment and installs the
Python modules directly into that environment's `site-packages`.

### 1. Activate the environment and install build tools

If you already have an OpenMM environment, activate it and add the build tools:

```bash
conda activate YOUR_OPENMM_ENV
conda install -c conda-forge cmake swig make cxx-compiler
```

For micromamba:

```bash
micromamba activate YOUR_OPENMM_ENV
micromamba install -c conda-forge cmake swig make cxx-compiler
```

Confirm that OpenMM is imported from the intended environment:

```bash
python -c "import openmm; print(openmm.__version__, openmm.__file__)"
```

### 2. Unpack the OpenMM and Colvars sources

Download source archives from:

- [OpenMM releases](https://github.com/openmm/openmm/releases)
- [Colvars repository](https://github.com/Colvars/colvars)

Create `tests_dependence_src` in this repository, extract both archives there,
and rename the extracted directories to `openmm` and `colvars`. The required
layout is:

```text
openmm-colvars/
├── CMakeLists.txt
├── scripts/
│   └── build_install_python_env.sh
└── tests_dependence_src/
    ├── openmm/
    │   ├── CMakeLists.txt
    │   ├── libraries/lepton/
    │   └── platforms/cuda/        # required for a CUDA build
    └── colvars/
        └── src/
            └── colvarmodule.cpp
```

The directory names matter because the script uses these paths by default. Do
not commit the extracted dependency trees; the root `.gitignore` excludes them.

### 3. Build and install

From the repository root, choose one platform configuration.

CPU/Reference build:

```bash
export OPENMM_COLVARS_ENABLE_CUDA_KERNELS=OFF
./scripts/build_install_python_env.sh
```

CUDA build:

```bash
export OPENMM_COLVARS_ENABLE_CUDA_KERNELS=ON
./scripts/build_install_python_env.sh
```

CUDA is enabled by default by the script, but setting the variable explicitly
makes the selected build reproducible. The script locates `OpenMMCUDA` under
`$CONDA_PREFIX` and looks for CUDA headers in the environment, `CUDA_HOME`,
`CUDA_PATH`, or `/usr/local/cuda`.

If automatic detection fails, specify the paths explicitly:

```bash
export OPENMM_CUDA_LIBRARY="$CONDA_PREFIX/lib/plugins/libOpenMMCUDA.so"
export OPENMM_COLVARS_CUDA_TOOLKIT_INCLUDE_DIR="/usr/local/cuda/include"
./scripts/build_install_python_env.sh
```

You may also override the default source and build directories:

```bash
export OPENMM_SOURCE_DIR="/path/to/openmm"
export COLVARS_SOURCE_DIR="/path/to/colvars"
export OPENMM_COLVARS_BUILD_DIR="/path/to/build"
./scripts/build_install_python_env.sh
```

The installer builds `OpenMMColvars`, generates the SWIG extension, and installs
`openmm_colvars`, `_openmm_colvars`, and `openmm_fep` into the active Python
environment. No `PYTHONPATH` modification is required.

## Quick verification

```bash
python - <<'PY'
from openmm import Platform, System
from openmm_colvars import ColvarsForce
from openmm_fep import FEPAtoms, FEPConfig

platforms = [
    Platform.getPlatform(i).getName()
    for i in range(Platform.getNumPlatforms())
]
print("OpenMM platforms:", platforms)

system = System()
system.addParticle(1.0)
system.addParticle(1.0)
force = ColvarsForce.fromConfigString(
    "colvarsTrajFrequency 0\n",
    "colvars_pycheck",
)
system.addForce(force)
print("Colvars force:", type(force).__name__)
print("FEP API:", FEPAtoms, FEPConfig)
PY
```

To confirm that CUDA itself is usable:

```bash
python - <<'PY'
from openmm import Platform
cuda = Platform.getPlatformByName("CUDA")
print("CUDA platform:", cuda.getName())
PY
```

## Manual CMake build

The installation script is recommended for Python users. The equivalent core
CMake configuration is:

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" \
  -DOPENMM_SOURCE_DIR="$PWD/tests_dependence_src/openmm" \
  -DCOLVARS_SOURCE_DIR="$PWD/tests_dependence_src/colvars" \
  -DLEPTON_DIR="$PWD/tests_dependence_src/openmm/libraries/lepton" \
  -DOPENMM_COLVARS_BUILD_PYTHON=ON \
  -DOPENMM_COLVARS_ENABLE_CUDA_KERNELS=OFF

cmake --build build -j
cmake --install build
```

Change `OPENMM_COLVARS_ENABLE_CUDA_KERNELS` to `ON` and provide the CUDA library
and include directory described above for a CUDA ABF/eABF build.

## Usage

### Minimal Python example

The following example creates a two-particle OpenMM system, adds a Colvars
distance restraint, and runs it on the CPU platform:

```python
from pathlib import Path

from openmm import (
    Context,
    HarmonicBondForce,
    Platform,
    System,
    Vec3,
    VerletIntegrator,
)
from openmm.unit import (
    dalton,
    kilojoule_per_mole,
    nanometer,
    picosecond,
)
from openmm_colvars import ColvarsForce

# Colvars writes trajectory/restart files using this prefix.
Path("output").mkdir(exist_ok=True)

system = System()
system.addParticle(39.9 * dalton)
system.addParticle(39.9 * dalton)

# Add physical forces before ColvarsForce.
bond = HarmonicBondForce()
bond.addBond(
    0,
    1,
    0.30 * nanometer,
    1000.0 * kilojoule_per_mole / nanometer**2,
)
system.addForce(bond)

# Colvars atomNumbers are one-based. Colvars length values use angstroms, so
# this 3.0 angstrom restraint corresponds to the 0.30 nm OpenMM bond length.
colvars_config = """
colvarsTrajFrequency 100

colvar {
  name d12
  width 0.1
  distance {
    group1 { atomNumbers 1 }
    group2 { atomNumbers 2 }
  }
}

harmonic {
  name restrain_d12
  colvars d12
  centers 3.0
  forceConstant 1.0
}
"""

colvars = ColvarsForce.fromConfigString(
    colvars_config,
    "output/colvars",
)
system.addForce(colvars)

integrator = VerletIntegrator(0.001 * picosecond)
platform = Platform.getPlatformByName("CPU")
context = Context(system, integrator, platform)
context.setPositions(
    [
        Vec3(0.0, 0.0, 0.0),
        Vec3(0.3, 0.0, 0.0),
    ]
    * nanometer
)

integrator.step(1000)
state = context.getState(getEnergy=True, getPositions=True)
print("Platform:", context.getPlatform().getName())
print("Potential energy:", state.getPotentialEnergy())
print("Positions:", state.getPositions())
```

To read the Colvars configuration from a file instead, replace
`fromConfigString(...)` with:

```python
colvars = ColvarsForce("colvars.in", "output/colvars")
```

Add `ColvarsForce` after the physical forces that should contribute to
ABF/eABF total-force collection.

### C++

```cpp
#include "ColvarsForce.h"

OpenMM::ColvarsForce colvars("colvars.in", "output/colvars");
colvars.setLambdaParameterName("lambda");
system.addForce(&colvars);
```

To construct a force from an in-memory configuration:

```cpp
auto colvars = OpenMM::ColvarsForce::fromConfigString(configText, "colvars");
```

### Alchemical lambda

`alchLambda` reads and writes an OpenMM context parameter. Its default name is
`lambda`; change it with `setLambdaParameterName()`.

`alchFlambda` requires `dE/dlambda`. The interface first reads OpenMM's internal
energy-parameter derivative cache. If no derivative is present, it falls back
to a context parameter named `dE_dlambda`; change that name with
`setLambdaDerivativeParameterName()`.

### TABF alpha

When the Colvars source provides `tabfAlpha`, the interface implements the TABF
backend hooks. Enable it on each force by calling `setTabfEnabled(true)` before
creating the OpenMM `Context`.

The default OpenMM variable names are `tabf_alpha` and `tabf_dE_dalpha`.
`ColvarsTABFIntegrator` applies physical forces as `f/alpha`, leaves Colvars
bias forces unscaled, and publishes
`d(U/alpha)/dalpha = -U/alpha^2`. Put physical forces and `ColvarsForce` in
separate force groups; the default groups are 0 and 1, respectively.

When combining `alchLambda` with TABF, call
`setScaleLambdaDerivativeByTabfAlpha(true)` so the proxy returns
`dU/dlambda / alpha`.

## OpenMM FEP helpers

The included `openmm_fep` package builds and runs dual-topology FEP windows with
NAMD-style lambda schedules, soft-core Lennard-Jones terms, and NAMD-like
`.fepout` output.

### Minimal example

```python
from openmm_fep import FEPAtoms, FEPConfig, prepare_fep_system, run_fep

config = FEPConfig(
    alch_atoms=FEPAtoms(vanish=[0, 1, 2], grow=[50, 51, 52]),
    # Or: FEPAtoms.from_alch_file("mutation.fep", alch_col="B")
    alch_vdw_shift_coeff=5.0,
    alch_elec_lambda_start=0.5,
    alch_vdw_lambda_end=1.0,
    windows=[0.0, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0],
    steps=[50000, 50000, 50000, 50000, 50000, 50000, 50000],
    alch_out_file="run.fepout",
    alch_out_freq=5,
    alch_equil_steps=0,
    alch_decouple=False,
    double_wide=True,
)

prepared = prepare_fep_system(system, config)

def simulation_factory(prepared_system):
    # Return an openmm.app.Simulation using prepared_system.
    return simulation

run_fep(prepared, simulation_factory)
```

### Implemented NAMD alchemy semantics

- `alchFile` is parsed as a PDB-like file. By default, beta values of `-1`, `1`,
  and `0` mean vanish, grow, and unchanged.
- Python atom indices are zero-based.
- Adjacent `windows` values define each `lambda/lambda2` window, so
  `len(steps) == len(windows) - 1`.
- `double_wide=True` adds NAMD's final single-direction endpoint window when a
  schedule ends at lambda 0 or 1.
- NAMD-style staggered electrostatic and van der Waals schedules are used.
- Soft-core van der Waals interactions use `alch_vdw_shift_coeff`, default 5.0.
- Grow-versus-vanish interactions are excluded for dual-topology end states.
- Alchemical 1-4 exceptions are handled separately for electrostatics and van
  der Waals interactions.
- Output uses NAMD-like headers, `FepEnergy:`, `FepE_back:`, and per-window
  free-energy footer lines.

## Units

Colvars/OpenMM uses the OpenMM unit system: nm, kJ/mol.

## Periodic boundary conditions

Periodic boundary conditions default to `PbcMode::Auto`: Colvars PBC is enabled
only when another force in the same OpenMM `System` declares periodic boundary
conditions. Use `PbcMode::Enabled` or `PbcMode::Disabled` to override this.

## Limitations

- Scripted Colvars functions and scripted forces are not supported.
- Multiple-replica communication is not supported.
- Volumetric map variables and GPU-resident Colvars kernels are not supported.
- OpenCL, HIP, and CUDA builds without the deferred CUDA path cannot perform
  same-step ABF/eABF total-force sampling.
- `openmm_fep` expects a standard OpenMM system with exactly one
  `NonbondedForce`.
- `alch_decouple=True` is not implemented.
- Arbitrary pre-existing `CustomNonbondedForce` objects require force-specific
  alchemical rules. OpenMM CHARMM NBFIX terms are handled explicitly.
- PME electrostatics use OpenMM-native charge scaling and are not intended to
  reproduce NAMD's internal PME decomposition byte-for-byte.
