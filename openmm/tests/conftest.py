"""Local OpenMM regression selection; no remote or CI execution is involved."""

from dataclasses import dataclass

import openmm as mm
import pytest


def pytest_addoption(parser):
    group = parser.getgroup("OpenMM-Colvars")
    group.addoption("--openmm-platform", default="Reference,CPU",
                    help="Comma-separated Reference,CPU,CUDA platforms")
    group.addoption("--openmm-precision", choices=("single", "mixed", "double"),
                    default="mixed", help="CUDA precision (Reference uses double)")
    group.addoption("--require-cuda", action="store_true",
                    help="Fail rather than skip if a requested CUDA device is unavailable")


@dataclass(frozen=True)
class Engine:
    platform: str
    use_gpu: bool
    precision: str

    @property
    def properties(self):
        if self.platform == "CUDA":
            return {"Precision": self.precision, "DeterministicForces": "true"}
        if self.platform == "CPU":
            return {"Threads": "1"}
        return {}

    @property
    def tolerances(self):
        # Absolute tolerances: CV nm, energy kJ/mol, force kJ/mol/nm.
        if self.platform == "CUDA" and self.precision == "single":
            return (2.0e-6, 2.0e-4, 2.0e-3)
        if self.platform == "CUDA" and self.precision == "mixed":
            return (2.0e-8, 2.0e-6, 2.0e-5)
        return (1.0e-12, 1.0e-9, 1.0e-8)


def pytest_generate_tests(metafunc):
    if "engine" not in metafunc.fixturenames:
        return
    params = []
    for name in metafunc.config.getoption("--openmm-platform").split(","):
        name = name.strip()
        if name not in ("Reference", "CPU", "CUDA"):
            raise pytest.UsageError(f"Unsupported test platform {name!r}")
        precision = metafunc.config.getoption("--openmm-precision") if name == "CUDA" else "double"
        modes = (True, False) if name == "CUDA" else (False,)
        for mode in modes:
            params.append(pytest.param(Engine(name, mode, precision),
                          id=f"{name}-{'gpu' if mode else 'cpu'}-{precision}"))
    metafunc.parametrize("engine", params, indirect=True)


@pytest.fixture
def engine(request):
    spec = request.param
    try:
        platform = mm.Platform.getPlatformByName(spec.platform)
    except mm.OpenMMException as exc:
        if spec.platform != "CUDA" or request.config.getoption("--require-cuda"):
            pytest.fail(f"Requested OpenMM platform {spec.platform} unavailable: {exc}")
        pytest.skip(f"CUDA platform unavailable: {exc}")
    if spec.platform == "CUDA":
        # Probe bare OpenMM separately, so a broken Colvars CUDA implementation
        # cannot be misreported as a missing GPU and silently skipped.
        system = mm.System()
        system.addParticle(1.0)
        integrator = mm.VerletIntegrator(0.001)
        try:
            context = mm.Context(system, integrator, platform, spec.properties)
        except mm.OpenMMException as exc:
            if request.config.getoption("--require-cuda"):
                pytest.fail(f"Bare OpenMM CUDA context unavailable: {exc}")
            pytest.skip(f"Bare OpenMM CUDA context unavailable: {exc}")
        else:
            del context
    return spec
