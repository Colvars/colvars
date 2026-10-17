from __future__ import annotations

import copy
from dataclasses import dataclass
import os
from pathlib import Path
import sys
from typing import Callable

from .config import FEPConfig, LambdaWindow
from .fepout import FEPOutWriter, FEPStatistics, kj_per_mol_to_kcal_per_mol

_OPENMM_PLUGINS_LOADED = False


@dataclass
class PreparedFEPSystem:
    system: object
    config: FEPConfig
    force_group_electrostatics: int
    force_group_vdw: int
    parameter_names: tuple[str, ...] = (
        "fep_grow_electrostatics",
        "fep_vanish_electrostatics",
        "fep_grow_electrostatics_sq",
        "fep_vanish_electrostatics_sq",
        "fep_grow_vdw",
        "fep_vanish_vdw",
    )


@dataclass
class PreparedLambdaABFSystem:
    system: object
    config: FEPConfig
    colvars_force: object | None
    force_group_electrostatics: int
    force_group_vdw: int
    lambda_parameter_name: str = "lambda"
    lambda_derivative_parameter_name: str = "dE_dlambda"
    alchemical_pme_force: object | None = None


def _require_openmm():
    try:
        import openmm  # type: ignore
    except ImportError as exc:
        raise ImportError("OpenMM is required for this operation") from exc
    return openmm


def load_openmm_platform_plugins() -> tuple[str, ...]:
    """Load OpenMM plugins from an unactivated conda/micromamba environment.

    On Windows, launching ``.../env/python.exe`` directly does not add
    ``Library/bin`` to the DLL search path.  CUDA/OpenCL plugins can then fail
    to register even though their DLLs exist in the environment.
    """

    global _OPENMM_PLUGINS_LOADED
    openmm = _require_openmm()
    if _OPENMM_PLUGINS_LOADED:
        return ()
    env_prefix = Path(sys.executable).resolve().parent
    candidates = (
        env_prefix / "Library" / "bin",
        env_prefix / "Library" / "lib",
        env_prefix / "Library" / "lib" / "plugins",
    )
    for path in candidates:
        if path.exists():
            if hasattr(os, "add_dll_directory"):
                os.add_dll_directory(str(path))
            os.environ["PATH"] = str(path) + os.pathsep + os.environ.get("PATH", "")

    plugin_dir = env_prefix / "Library" / "lib" / "plugins"
    if plugin_dir.exists():
        existing = {openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())}
        cuda_plugin = plugin_dir / "OpenMMCUDA.dll"
        if "CUDA" not in existing and cuda_plugin.exists():
            try:
                openmm.Platform.loadPluginLibrary(str(cuda_plugin))
            except Exception:
                pass
            _OPENMM_PLUGINS_LOADED = True
            return ()
        try:
            loaded = tuple(openmm.Platform.loadPluginsFromDirectory(str(plugin_dir)))
            _OPENMM_PLUGINS_LOADED = True
            return loaded
        except Exception:
            _OPENMM_PLUGINS_LOADED = True
            return ()
    _OPENMM_PLUGINS_LOADED = True
    return ()


def available_platform_names() -> tuple[str, ...]:
    openmm = _require_openmm()
    load_openmm_platform_plugins()
    return tuple(openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms()))


def select_openmm_platform(device_index: str | int | None = None):
    openmm = _require_openmm()
    names = available_platform_names()
    if "CUDA" in names:
        platform = openmm.Platform.getPlatformByName("CUDA")
        properties = {}
        if device_index is not None:
            properties["DeviceIndex"] = str(device_index)
        return platform, properties
    return openmm.Platform.getPlatformByName("Reference"), {}


def _quantity_value(value, unit_value) -> float:
    if hasattr(value, "value_in_unit"):
        return float(value.value_in_unit(unit_value))
    return float(value)


def _charge_value(value) -> float:
    openmm = _require_openmm()
    return _quantity_value(value, openmm.unit.elementary_charge)


def _charge_product_value(value) -> float:
    openmm = _require_openmm()
    return _quantity_value(value, openmm.unit.elementary_charge**2)


def _sigma_value(value) -> float:
    openmm = _require_openmm()
    return _quantity_value(value, openmm.unit.nanometer)


def _epsilon_value(value) -> float:
    openmm = _require_openmm()
    return _quantity_value(value, openmm.unit.kilojoule_per_mole)


def _copy_force_settings(source, target) -> None:
    openmm = _require_openmm()
    method = source.getNonbondedMethod()
    if method == openmm.NonbondedForce.NoCutoff:
        target.setNonbondedMethod(openmm.CustomNonbondedForce.NoCutoff)
    elif method == openmm.NonbondedForce.CutoffNonPeriodic:
        target.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffNonPeriodic)
    else:
        target.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
    if hasattr(source, "getCutoffDistance"):
        target.setCutoffDistance(source.getCutoffDistance())
    if hasattr(source, "getUseSwitchingFunction"):
        target.setUseSwitchingFunction(source.getUseSwitchingFunction())
        target.setSwitchingDistance(source.getSwitchingDistance())
    if hasattr(source, "getUseLongRangeCorrection"):
        target.setUseLongRangeCorrection(source.getUseLongRangeCorrection())
    elif hasattr(source, "getUseDispersionCorrection"):
        target.setUseLongRangeCorrection(source.getUseDispersionCorrection())


def _copy_pme_force_settings(source, target) -> None:
    openmm = _require_openmm()
    method = source.getNonbondedMethod()
    if method == openmm.NonbondedForce.NoCutoff:
        target.setNonbondedMethod(target.NoCutoff)
    elif method in (openmm.NonbondedForce.CutoffPeriodic, openmm.NonbondedForce.CutoffNonPeriodic):
        if method == openmm.NonbondedForce.CutoffNonPeriodic:
            raise NotImplementedError("lambda-ABF alchemical PME does not support CutoffNonPeriodic")
        target.setNonbondedMethod(target.CutoffPeriodic)
    elif method == openmm.NonbondedForce.PME:
        target.setNonbondedMethod(target.PME)
        alpha, nx, ny, nz = source.getPMEParameters()
        alpha_nm = _quantity_value(alpha, openmm.unit.nanometer**-1)
        if alpha_nm <= 0.0 or nx <= 0 or ny <= 0 or nz <= 0:
            raise ValueError(
                "lambda-ABF alchemical PME requires explicit PME parameters; "
                "call NonbondedForce.setPMEParameters(alpha, nx, ny, nz) first"
            )
        target.setPMEParameters(alpha_nm, int(nx), int(ny), int(nz))
    else:
        raise NotImplementedError("lambda-ABF alchemical PME requires NoCutoff, CutoffPeriodic, or PME")
    if hasattr(source, "getCutoffDistance"):
        cutoff_nm = _quantity_value(source.getCutoffDistance(), openmm.unit.nanometer)
        target.setCutoffDistance(cutoff_nm)
    if hasattr(source, "getReactionFieldDielectric") and hasattr(target, "setReactionFieldDielectric"):
        target.setReactionFieldDielectric(source.getReactionFieldDielectric())


def _force_uses_periodic_nonbonded(force) -> bool:
    openmm = _require_openmm()
    if isinstance(force, openmm.NonbondedForce):
        return force.getNonbondedMethod() in (
            openmm.NonbondedForce.CutoffPeriodic,
            openmm.NonbondedForce.PME,
            openmm.NonbondedForce.LJPME,
        )
    if isinstance(force, openmm.CustomNonbondedForce):
        return force.getNonbondedMethod() == openmm.CustomNonbondedForce.CutoffPeriodic
    if hasattr(force, "usesPeriodicBoundaryConditions"):
        return bool(force.usesPeriodicBoundaryConditions())
    return False


def _system_uses_periodic_nonbonded(system) -> bool:
    return any(_force_uses_periodic_nonbonded(system.getForce(i)) for i in range(system.getNumForces()))


def _find_free_force_group(system, reserved: set[int] | None = None) -> int:
    reserved = reserved or set()
    used = {system.getForce(i).getForceGroup() for i in range(system.getNumForces())}
    for group in range(32):
        if group not in used and group not in reserved:
            return group
    raise ValueError("no free OpenMM force groups available")


def _is_nonbonded_exception(nbforce, atom_i: int, atom_j: int) -> bool:
    wanted = {atom_i, atom_j}
    for idx in range(nbforce.getNumExceptions()):
        i, j, *_ = nbforce.getExceptionParameters(idx)
        if {int(i), int(j)} == wanted:
            return True
    return False


def _exclusion_key(atom_i: int, atom_j: int) -> tuple[int, int]:
    return (atom_i, atom_j) if atom_i < atom_j else (atom_j, atom_i)


def _has_custom_exclusion(force, atom_i: int, atom_j: int) -> bool:
    wanted = _exclusion_key(atom_i, atom_j)
    for idx in range(force.getNumExclusions()):
        i, j = force.getExclusionParticles(idx)
        if _exclusion_key(int(i), int(j)) == wanted:
            return True
    return False


def _add_custom_exclusion_once(force, atom_i: int, atom_j: int) -> None:
    if not _has_custom_exclusion(force, atom_i, atom_j):
        force.addExclusion(atom_i, atom_j)


def _has_nonexcluded_pair_within(nbforce, atoms: set[int]) -> bool:
    if len(atoms) < 2:
        return False
    exclusions = set()
    for idx in range(nbforce.getNumExceptions()):
        i, j, *_ = nbforce.getExceptionParameters(idx)
        exclusions.add(_exclusion_key(int(i), int(j)))
    ordered = sorted(atoms)
    for offset, atom_i in enumerate(ordered):
        for atom_j in ordered[offset + 1 :]:
            if _exclusion_key(atom_i, atom_j) not in exclusions:
                return True
    return False


def _has_nonexcluded_pair_between(nbforce, atoms1: set[int], atoms2: set[int]) -> bool:
    if not atoms1 or not atoms2:
        return False
    exclusions = set()
    for idx in range(nbforce.getNumExceptions()):
        i, j, *_ = nbforce.getExceptionParameters(idx)
        exclusions.add(_exclusion_key(int(i), int(j)))
    for atom_i in atoms1:
        for atom_j in atoms2:
            if atom_i == atom_j:
                continue
            if _exclusion_key(atom_i, atom_j) not in exclusions:
                return True
    return False


def _is_charmm_nbfix_force(force) -> bool:
    expression = force.getEnergyFunction().replace(" ", "")
    return (
        "(a/r6)^2-b/r6" in expression
        and "a=acoef(type1,type2)" in expression
        and "b=bcoef(type1,type2)" in expression
        and force.getNumPerParticleParameters() == 1
        and force.getPerParticleParameterName(0) == "type"
        and {force.getTabulatedFunctionName(i) for i in range(force.getNumTabulatedFunctions())}
        >= {"acoef", "bcoef"}
    )


def _alchemically_modify_nbfix_force(force, config: FEPConfig, grow: set[int], vanish: set[int], group_vdw: int) -> None:
    shift_nm2 = float(config.alch_vdw_shift_coeff) * 0.01
    grow_pair = "step(is_grow1+is_grow2-0.5)"
    vanish_pair = "step(is_vanish1+is_vanish2-0.5)"
    mixed_pair = f"({grow_pair})*({vanish_pair})"
    lam_pair = (
        f"(1-({mixed_pair}))*select(({grow_pair}), fep_grow_vdw, "
        f"select(({vanish_pair}), fep_vanish_vdw, 1))"
    )
    a_pair = "acoef(type1,type2)"
    b_pair = "bcoef(type1,type2)"
    rsoft = f"sqrt(r^2+{shift_nm2}*(1-({lam_pair})))"
    expression = f"({lam_pair})*((({a_pair})/(({rsoft})^6))^2-({b_pair})/(({rsoft})^6))"
    force.addGlobalParameter("fep_grow_vdw", 1.0)
    force.addGlobalParameter("fep_vanish_vdw", 1.0)
    force.addPerParticleParameter("is_grow")
    force.addPerParticleParameter("is_vanish")
    for atom in range(force.getNumParticles()):
        params = list(force.getParticleParameters(atom))
        params.extend([1.0 if atom in grow else 0.0, 1.0 if atom in vanish else 0.0])
        force.setParticleParameters(atom, params)
    force.setEnergyFunction(expression)
    force.setForceGroup(group_vdw)


def _partition(atom: int, grow: set[int], vanish: set[int]) -> str:
    if atom in grow:
        return "grow"
    if atom in vanish:
        return "vanish"
    return "env"


def _exception_parameter_name(atom_i: int, atom_j: int, grow: set[int], vanish: set[int]) -> str | None:
    pi = _partition(atom_i, grow, vanish)
    pj = _partition(atom_j, grow, vanish)
    parts = {pi, pj}
    if parts == {"grow", "vanish"}:
        return None
    if parts == {"grow", "env"}:
        return "fep_grow_electrostatics"
    if parts == {"vanish", "env"}:
        return "fep_vanish_electrostatics"
    if parts == {"grow"}:
        return "fep_grow_electrostatics"
    if parts == {"vanish"}:
        return "fep_vanish_electrostatics"
    return None


def _add_same_partition_electrostatic_correction(
    system, nbforce, original_particle_params, config: FEPConfig, grow: set[int], vanish: set[int], group_elec: int
) -> None:
    openmm = _require_openmm()
    add_grow_group = _has_nonexcluded_pair_within(nbforce, grow)
    add_vanish_group = _has_nonexcluded_pair_within(nbforce, vanish)
    if config.alch_decouple or not (add_grow_group or add_vanish_group):
        return

    expression = (
        "138.935456*charge1*charge2*lambda_correction/r;"
        "lambda_correction=select(is_grow1+is_grow2, "
        "fep_grow_electrostatics-fep_grow_electrostatics_sq, "
        "fep_vanish_electrostatics-fep_vanish_electrostatics_sq)"
    )
    correction = openmm.CustomNonbondedForce(expression)
    correction.addGlobalParameter("fep_grow_electrostatics", 1.0)
    correction.addGlobalParameter("fep_vanish_electrostatics", 1.0)
    correction.addGlobalParameter("fep_grow_electrostatics_sq", 1.0)
    correction.addGlobalParameter("fep_vanish_electrostatics_sq", 1.0)
    correction.addPerParticleParameter("charge")
    correction.addPerParticleParameter("is_grow")
    for atom in range(system.getNumParticles()):
        charge, _sigma, _epsilon = original_particle_params[atom]
        correction.addParticle([charge, 1.0 if atom in grow else 0.0])
    _copy_force_settings(nbforce, correction)
    correction.setForceGroup(group_elec)
    if add_grow_group:
        correction.addInteractionGroup(grow, grow)
    if add_vanish_group:
        correction.addInteractionGroup(vanish, vanish)
    for idx in range(nbforce.getNumExceptions()):
        i, j, *_ = nbforce.getExceptionParameters(idx)
        correction.addExclusion(int(i), int(j))
    system.addForce(correction)


def _lambda_schedule_expression(
    name: str, config: FEPConfig, lambda_parameter_name: str, clamp_lambda: bool = True
) -> str:
    lam = f"min(1,max(0,{lambda_parameter_name}))" if clamp_lambda else lambda_parameter_name
    if name == "elec_grow":
        start = float(config.alch_elec_lambda_start)
        if start >= 1.0:
            return f"step(({lam})-1)"
        return f"select(step(({lam})-{start}), (({lam})-{start})/{1.0-start}, 0)"
    if name == "elec_vanish":
        start = float(config.alch_elec_lambda_start)
        if start >= 1.0:
            return f"1-step({lam})"
        end = 1.0 - start
        return f"select(step({end}-({lam})), 1-({lam})/{end}, 0)"
    if name == "vdw_grow":
        end = float(config.alch_vdw_lambda_end)
        if end <= 0.0:
            return "1"
        return f"select(step({end}-({lam})), ({lam})/{end}, 1)"
    if name == "vdw_vanish":
        end = float(config.alch_vdw_lambda_end)
        if end <= 0.0:
            return "0"
        start = 1.0 - end
        return f"select(step(({lam})-{start}), (1-({lam}))/{end}, 1)"
    raise ValueError(f"unknown lambda schedule {name!r}")


def _add_same_partition_electrostatic_lambda_correction(
    system,
    nbforce,
    original_particle_params,
    config: FEPConfig,
    grow: set[int],
    vanish: set[int],
    group_elec: int,
    lambda_parameter_name: str,
) -> None:
    openmm = _require_openmm()
    add_grow_group = _has_nonexcluded_pair_within(nbforce, grow)
    add_vanish_group = _has_nonexcluded_pair_within(nbforce, vanish)
    if config.alch_decouple or not (add_grow_group or add_vanish_group):
        return

    lam = "lambda_clamped"
    grow_elec = _lambda_schedule_expression("elec_grow", config, lam, clamp_lambda=False)
    vanish_elec = _lambda_schedule_expression("elec_vanish", config, lam, clamp_lambda=False)
    expression = (
        "138.935456*charge1*charge2*lambda_correction/r;"
        "lambda_correction=select(is_grow1+is_grow2, "
        "grow_elec-grow_elec^2, vanish_elec-vanish_elec^2);"
        f"grow_elec={grow_elec};"
        f"vanish_elec={vanish_elec};"
        f"lambda_clamped=min(1,max(0,{lambda_parameter_name}))"
    )
    correction = openmm.CustomNonbondedForce(expression)
    correction.addGlobalParameter(lambda_parameter_name, 0.0)
    correction.addEnergyParameterDerivative(lambda_parameter_name)
    correction.addPerParticleParameter("charge")
    correction.addPerParticleParameter("is_grow")
    for atom in range(system.getNumParticles()):
        charge, _sigma, _epsilon = original_particle_params[atom]
        correction.addParticle([charge, 1.0 if atom in grow else 0.0])
    _copy_force_settings(nbforce, correction)
    correction.setForceGroup(group_elec)
    if add_grow_group:
        correction.addInteractionGroup(grow, grow)
    if add_vanish_group:
        correction.addInteractionGroup(vanish, vanish)
    for idx in range(nbforce.getNumExceptions()):
        i, j, *_ = nbforce.getExceptionParameters(idx)
        correction.addExclusion(int(i), int(j))
    system.addForce(correction)


def _alchemical_pme_charge_endpoints(atom: int, charge: float, grow: set[int], vanish: set[int]) -> tuple[float, float]:
    if atom in grow:
        return 0.0, charge
    if atom in vanish:
        return charge, 0.0
    return charge, charge


def _alchemical_pme_exception_endpoints(
    atom_i: int,
    atom_j: int,
    charge_prod: float,
    grow: set[int],
    vanish: set[int],
) -> tuple[float, float]:
    pi = _partition(atom_i, grow, vanish)
    pj = _partition(atom_j, grow, vanish)
    parts = {pi, pj}
    if parts == {"grow", "vanish"}:
        return 0.0, 0.0
    if "grow" in parts and parts <= {"grow", "env"}:
        return 0.0, charge_prod
    if "vanish" in parts and parts <= {"vanish", "env"}:
        return charge_prod, 0.0
    return charge_prod, charge_prod


def prepare_lambda_abf_system(
    system,
    config: FEPConfig,
    colvars_force=None,
    lambda_parameter_name: str = "lambda",
    lambda_derivative_parameter_name: str = "dE_dlambda",
) -> PreparedLambdaABFSystem:
    """Prepare a continuous-lambda OpenMM System for Colvars alchLambda/eABF.

    The returned System is a copy.  Electrostatics are moved into
    ``openmm_colvars.AlchemicalPMEForce`` so PME energy, forces and
    ``dE/dlambda`` are available to Colvars.  Soft-core vdW terms are custom
    forces whose expressions depend directly on the same OpenMM global
    ``lambda`` parameter.
    """

    openmm = _require_openmm()
    try:
        from openmm_colvars import AlchemicalPMEForce
    except ImportError as exc:
        raise ImportError(
            "prepare_lambda_abf_system requires the openmm_colvars Python extension"
        ) from exc

    if config.alch_decouple:
        raise NotImplementedError("alch_decouple=True is not implemented yet")
    new_system = copy.deepcopy(system)
    config.validate(new_system.getNumParticles())

    nbforces = []
    custom_nonbonded = []
    for i in range(new_system.getNumForces()):
        force = new_system.getForce(i)
        if isinstance(force, openmm.NonbondedForce):
            nbforces.append(force)
        elif isinstance(force, openmm.CustomNonbondedForce):
            custom_nonbonded.append(force)
    if len(nbforces) != 1:
        raise ValueError("prepare_lambda_abf_system requires exactly one OpenMM NonbondedForce")
    if custom_nonbonded:
        raise NotImplementedError(
            "prepare_lambda_abf_system does not support existing CustomNonbondedForce terms yet"
        )

    nbforce = nbforces[0]
    group_elec = _find_free_force_group(new_system)
    group_vdw = _find_free_force_group(new_system, {group_elec})
    grow = set(config.alch_atoms.grow)
    vanish = set(config.alch_atoms.vanish)
    alchemical = grow | vanish
    chemical = set(range(new_system.getNumParticles())) - alchemical
    original_particle_params = {
        atom: nbforce.getParticleParameters(atom) for atom in range(new_system.getNumParticles())
    }

    alch_pme = AlchemicalPMEForce()
    alch_pme.setLambdaParameterName(lambda_parameter_name)
    alch_pme.setLambdaDerivativeParameterName(lambda_derivative_parameter_name)
    if hasattr(alch_pme, "setElectrostaticLambdaStart"):
        alch_pme.setElectrostaticLambdaStart(float(config.alch_elec_lambda_start))
    _copy_pme_force_settings(nbforce, alch_pme)
    alch_pme.setForceGroup(group_elec)

    for atom in range(new_system.getNumParticles()):
        charge, sigma, epsilon = original_particle_params[atom]
        charge_float = _charge_value(charge)
        sigma_float = _sigma_value(sigma)
        charge0, charge1 = _alchemical_pme_charge_endpoints(atom, charge_float, grow, vanish)
        alch_pme.addParticle(charge0, charge1, sigma_float, sigma_float, 0.0, 0.0)

        if atom in alchemical:
            nbforce.setParticleParameters(atom, 0.0 * charge, sigma, 0.0 * epsilon)
        else:
            nbforce.setParticleParameters(atom, 0.0 * charge, sigma, epsilon)

    exception_lj_bonds: list[tuple[int, int, object, object]] = []
    for idx in range(nbforce.getNumExceptions()):
        i, j, charge_prod, sigma, epsilon = nbforce.getExceptionParameters(idx)
        i = int(i)
        j = int(j)
        charge_prod_float = _charge_product_value(charge_prod)
        sigma_float = _sigma_value(sigma)
        charge_prod0, charge_prod1 = _alchemical_pme_exception_endpoints(
            i, j, charge_prod_float, grow, vanish
        )
        alch_pme.addException(i, j, charge_prod0, charge_prod1, sigma_float, sigma_float, 0.0, 0.0)

        if i in alchemical or j in alchemical:
            nbforce.setExceptionParameters(idx, i, j, 0.0 * charge_prod, sigma, 0.0 * epsilon)
            if not (i in grow and j in vanish) and not (i in vanish and j in grow):
                exception_lj_bonds.append((i, j, sigma, epsilon))
        else:
            nbforce.setExceptionParameters(idx, i, j, 0.0 * charge_prod, sigma, epsilon)

    for i in grow:
        for j in vanish:
            if not _is_nonbonded_exception(nbforce, i, j):
                nbforce.addException(i, j, 0.0, 1.0, 0.0, True)
                alch_pme.addException(i, j, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0)

    nbforce.setForceGroup(group_vdw)
    new_system.addForce(alch_pme)

    _add_same_partition_electrostatic_lambda_correction(
        new_system,
        nbforce,
        original_particle_params,
        config,
        grow,
        vanish,
        group_elec,
        lambda_parameter_name,
    )

    add_env_group = _has_nonexcluded_pair_between(nbforce, alchemical, chemical)
    add_grow_group = _has_nonexcluded_pair_within(nbforce, grow)
    add_vanish_group = _has_nonexcluded_pair_within(nbforce, vanish)
    if add_env_group or add_grow_group or add_vanish_group:
        shift_nm2 = float(config.alch_vdw_shift_coeff) * 0.01
        lambda_grow = _lambda_schedule_expression("vdw_grow", config, "lambda_clamped", clamp_lambda=False)
        lambda_vanish = _lambda_schedule_expression("vdw_vanish", config, "lambda_clamped", clamp_lambda=False)
        expression = (
            "lambda_pair*4*epsilon_pair*((sigma_pair/reff_pair)^12-(sigma_pair/reff_pair)^6);"
            f"reff_pair=sqrt(r^2+{shift_nm2}*(1-lambda_pair));"
            "lambda_pair=select(is_grow1+is_grow2, lambda_grow, lambda_vanish);"
            "sigma_pair=0.5*(sigma1+sigma2);"
            "epsilon_pair=sqrt(epsilon1*epsilon2);"
            f"lambda_grow={lambda_grow};"
            f"lambda_vanish={lambda_vanish};"
            f"lambda_clamped=min(1,max(0,{lambda_parameter_name}))"
        )
        custom = openmm.CustomNonbondedForce(expression)
        custom.addGlobalParameter(lambda_parameter_name, 0.0)
        custom.addEnergyParameterDerivative(lambda_parameter_name)
        custom.addPerParticleParameter("sigma")
        custom.addPerParticleParameter("epsilon")
        custom.addPerParticleParameter("is_grow")
        for atom in range(new_system.getNumParticles()):
            _charge, sigma, epsilon = original_particle_params[atom]
            custom.addParticle([sigma, epsilon, 1.0 if atom in grow else 0.0])
        _copy_force_settings(nbforce, custom)
        custom.setForceGroup(group_vdw)
        if add_env_group:
            custom.addInteractionGroup(alchemical, chemical)
        if not config.alch_decouple:
            if add_grow_group:
                custom.addInteractionGroup(grow, grow)
            if add_vanish_group:
                custom.addInteractionGroup(vanish, vanish)
        for idx in range(nbforce.getNumExceptions()):
            i, j, *_ = nbforce.getExceptionParameters(idx)
            custom.addExclusion(int(i), int(j))
        new_system.addForce(custom)

    if exception_lj_bonds:
        shift_nm2 = float(config.alch_vdw_shift_coeff) * 0.01
        lambda_grow = _lambda_schedule_expression("vdw_grow", config, "lambda_clamped", clamp_lambda=False)
        lambda_vanish = _lambda_schedule_expression("vdw_vanish", config, "lambda_clamped", clamp_lambda=False)
        bond_expression = (
            "lambda_bond*4*epsilon*((sigma/reff_bond)^12-(sigma/reff_bond)^6);"
            f"reff_bond=sqrt(r^2+{shift_nm2}*(1-lambda_bond));"
            "lambda_bond=lambda_type*lambda_grow+(1-lambda_type)*lambda_vanish;"
            f"lambda_grow={lambda_grow};"
            f"lambda_vanish={lambda_vanish};"
            f"lambda_clamped=min(1,max(0,{lambda_parameter_name}))"
        )
        bonds = openmm.CustomBondForce(bond_expression)
        bonds.addGlobalParameter(lambda_parameter_name, 0.0)
        bonds.addEnergyParameterDerivative(lambda_parameter_name)
        bonds.addPerBondParameter("sigma")
        bonds.addPerBondParameter("epsilon")
        bonds.addPerBondParameter("lambda_type")
        for i, j, sigma, epsilon in exception_lj_bonds:
            pname = _exception_parameter_name(i, j, grow, vanish)
            selector = 1.0 if pname and "grow" in pname else 0.0
            bonds.addBond(i, j, [sigma, epsilon, selector])
        bonds.setUsesPeriodicBoundaryConditions(_system_uses_periodic_nonbonded(new_system))
        bonds.setForceGroup(group_vdw)
        new_system.addForce(bonds)

    if colvars_force is not None:
        if hasattr(colvars_force, "setLambdaParameterName"):
            colvars_force.setLambdaParameterName(lambda_parameter_name)
        if hasattr(colvars_force, "setLambdaDerivativeParameterName"):
            colvars_force.setLambdaDerivativeParameterName(lambda_derivative_parameter_name)
        if hasattr(colvars_force, "setAlchemyEnabled"):
            colvars_force.setAlchemyEnabled(True)
        if hasattr(colvars_force, "setUsesPeriodicBoundaryConditions"):
            colvars_force.setUsesPeriodicBoundaryConditions(_system_uses_periodic_nonbonded(new_system))
        new_system.addForce(colvars_force)

    return PreparedLambdaABFSystem(
        new_system,
        config,
        colvars_force,
        group_elec,
        group_vdw,
        lambda_parameter_name,
        lambda_derivative_parameter_name,
        alch_pme,
    )


def prepare_fep_system(system, config: FEPConfig) -> PreparedFEPSystem:
    """Return a copied OpenMM System with NAMD-style alchemical nonbonded forces.

    Standard ``NonbondedForce`` systems are supported. OpenMM CHARMM NBFIX
    ``CustomNonbondedForce`` terms are recognized and alchemically rewritten;
    unrecognized custom nonbonded terms are rejected rather than guessed.
    """

    openmm = _require_openmm()
    if config.alch_decouple:
        raise NotImplementedError("alch_decouple=True is not implemented yet")
    new_system = copy.deepcopy(system)
    config.validate(new_system.getNumParticles())
    nbforces = []
    nbfix_forces = []
    unsupported_custom_nonbonded = []
    for i in range(new_system.getNumForces()):
        force = new_system.getForce(i)
        if isinstance(force, openmm.NonbondedForce):
            nbforces.append(force)
        elif isinstance(force, openmm.CustomNonbondedForce):
            if _is_charmm_nbfix_force(force):
                nbfix_forces.append(force)
            else:
                unsupported_custom_nonbonded.append(force)
    if len(nbforces) != 1:
        raise ValueError("prepare_fep_system requires exactly one OpenMM NonbondedForce")
    if unsupported_custom_nonbonded:
        raise NotImplementedError(
            "prepare_fep_system does not support systems that already contain "
            "unrecognized CustomNonbondedForce terms; recognized OpenMM CHARMM NBFIX "
            "forces are supported"
        )
    nbforce = nbforces[0]

    group_elec = _find_free_force_group(new_system)
    group_vdw = _find_free_force_group(new_system, {group_elec})
    nbforce.setForceGroup(group_elec)

    grow = set(config.alch_atoms.grow)
    vanish = set(config.alch_atoms.vanish)
    alchemical = grow | vanish
    chemical = set(range(new_system.getNumParticles())) - alchemical
    original_particle_params = {
        atom: nbforce.getParticleParameters(atom) for atom in range(new_system.getNumParticles())
    }

    for name in (
        "fep_grow_electrostatics",
        "fep_vanish_electrostatics",
        "fep_grow_electrostatics_sq",
        "fep_vanish_electrostatics_sq",
    ):
        nbforce.addGlobalParameter(name, 1.0)

    for atom in grow:
        charge, sigma, epsilon = nbforce.getParticleParameters(atom)
        nbforce.setParticleParameters(atom, 0.0 * charge, sigma, 0.0 * epsilon)
        nbforce.addParticleParameterOffset("fep_grow_electrostatics", atom, charge, 0.0 * sigma, 0.0 * epsilon)
    for atom in vanish:
        charge, sigma, epsilon = nbforce.getParticleParameters(atom)
        nbforce.setParticleParameters(atom, 0.0 * charge, sigma, 0.0 * epsilon)
        nbforce.addParticleParameterOffset("fep_vanish_electrostatics", atom, charge, 0.0 * sigma, 0.0 * epsilon)

    exception_lj_bonds: list[tuple[int, int, object, object]] = []
    original_exception_count = nbforce.getNumExceptions()
    for idx in range(original_exception_count):
        i, j, charge_prod, sigma, epsilon = nbforce.getExceptionParameters(idx)
        i = int(i)
        j = int(j)
        if i not in alchemical and j not in alchemical:
            continue
        if {i, j}.intersection(grow) and {i, j}.intersection(vanish):
            nbforce.setExceptionParameters(idx, i, j, 0.0 * charge_prod, sigma, 0.0 * epsilon)
            continue
        pname = _exception_parameter_name(i, j, grow, vanish)
        nbforce.setExceptionParameters(idx, i, j, 0.0 * charge_prod, sigma, 0.0 * epsilon)
        if pname is not None:
            nbforce.addExceptionParameterOffset(pname, idx, charge_prod, 0.0 * sigma, 0.0 * epsilon)
            exception_lj_bonds.append((i, j, sigma, epsilon))

    # Make grow and vanish atoms invisible to one another in the original force, even if
    # the input topology did not exclude the dual-topology cross pairs.
    for i in grow:
        for j in vanish:
            if not _is_nonbonded_exception(nbforce, i, j):
                nbforce.addException(i, j, 0.0, 1.0, 0.0, True)
            for nbfix_force in nbfix_forces:
                _add_custom_exclusion_once(nbfix_force, i, j)

    _add_same_partition_electrostatic_correction(
        new_system, nbforce, original_particle_params, config, grow, vanish, group_elec
    )

    if nbfix_forces:
        for nbfix_force in nbfix_forces:
            _alchemically_modify_nbfix_force(nbfix_force, config, grow, vanish, group_vdw)
    else:
        add_env_group = _has_nonexcluded_pair_between(nbforce, alchemical, chemical)
        add_grow_group = _has_nonexcluded_pair_within(nbforce, grow)
        add_vanish_group = _has_nonexcluded_pair_within(nbforce, vanish)
        if add_env_group or add_grow_group or add_vanish_group:
            shift_nm2 = float(config.alch_vdw_shift_coeff) * 0.01
            lam_pair = "select(is_grow1+is_grow2, fep_grow_vdw, fep_vanish_vdw)"
            sigma_pair = "0.5*(sigma1+sigma2)"
            epsilon_pair = "sqrt(epsilon1*epsilon2)"
            reff_pair = f"sqrt(r^2+{shift_nm2}*(1-({lam_pair})))"
            expression = (
                f"({lam_pair})*4*({epsilon_pair})*"
                f"((({sigma_pair})/({reff_pair}))^12-(({sigma_pair})/({reff_pair}))^6)"
            )
            custom = openmm.CustomNonbondedForce(expression)
            custom.addGlobalParameter("fep_grow_vdw", 1.0)
            custom.addGlobalParameter("fep_vanish_vdw", 1.0)
            custom.addPerParticleParameter("sigma")
            custom.addPerParticleParameter("epsilon")
            custom.addPerParticleParameter("is_grow")
            for atom in range(new_system.getNumParticles()):
                _charge, sigma, epsilon = original_particle_params[atom]
                custom.addParticle([sigma, epsilon, 1.0 if atom in grow else 0.0])
            _copy_force_settings(nbforce, custom)
            custom.setForceGroup(group_vdw)
            if add_env_group:
                custom.addInteractionGroup(alchemical, chemical)
            if not config.alch_decouple:
                if add_grow_group:
                    custom.addInteractionGroup(grow, grow)
                if add_vanish_group:
                    custom.addInteractionGroup(vanish, vanish)
            for idx in range(nbforce.getNumExceptions()):
                i, j, *_ = nbforce.getExceptionParameters(idx)
                custom.addExclusion(int(i), int(j))
            new_system.addForce(custom)

    if exception_lj_bonds:
        shift_nm2 = float(config.alch_vdw_shift_coeff) * 0.01
        lam_bond = "(lambda_type*fep_grow_vdw+(1-lambda_type)*fep_vanish_vdw)"
        reff_bond = f"sqrt(r^2+{shift_nm2}*(1-{lam_bond}))"
        bond_expression = (
            f"{lam_bond}*4*epsilon*((sigma/({reff_bond}))^12-(sigma/({reff_bond}))^6)"
        )
        bonds = openmm.CustomBondForce(bond_expression)
        bonds.addGlobalParameter("fep_grow_vdw", 1.0)
        bonds.addGlobalParameter("fep_vanish_vdw", 1.0)
        bonds.addPerBondParameter("sigma")
        bonds.addPerBondParameter("epsilon")
        bonds.addPerBondParameter("lambda_type")
        for i, j, sigma, epsilon in exception_lj_bonds:
            pname = _exception_parameter_name(i, j, grow, vanish)
            lambda_type = "fep_grow_vdw" if pname and "grow" in pname else "fep_vanish_vdw"
            # CustomBondForce per-bond parameters cannot directly select a global by name, so
            # split grow and vanish values by encoding a selector and using both globals.
            selector = 1.0 if lambda_type == "fep_grow_vdw" else 0.0
            bonds.addBond(i, j, [sigma, epsilon, selector])
        bonds.setUsesPeriodicBoundaryConditions(_system_uses_periodic_nonbonded(new_system))
        bonds.setForceGroup(group_vdw)
        new_system.addForce(bonds)

    return PreparedFEPSystem(new_system, config, group_elec, group_vdw)


def set_context_lambda(context, config: FEPConfig, lam: float) -> None:
    state = config.lambda_state(lam)
    available = context.getParameters()
    parameter_values = {
        "fep_grow_electrostatics": state.grow_electrostatics,
        "fep_vanish_electrostatics": state.vanish_electrostatics,
        "fep_grow_electrostatics_sq": state.grow_electrostatics * state.grow_electrostatics,
        "fep_vanish_electrostatics_sq": state.vanish_electrostatics * state.vanish_electrostatics,
        "fep_grow_vdw": state.grow_vdw,
        "fep_vanish_vdw": state.vanish_vdw,
    }
    for name, value in parameter_values.items():
        if name in available:
            context.setParameter(name, value)


def _energy_kj(context, groups: set[int] | None = None) -> float:
    state = context.getState(getEnergy=True, groups=groups) if groups is not None else context.getState(getEnergy=True)
    return float(state.getPotentialEnergy().value_in_unit_system(_require_openmm().unit.md_unit_system))


def _component_energies_kcal(context, prepared: PreparedFEPSystem) -> tuple[float, float]:
    elec = kj_per_mol_to_kcal_per_mol(_energy_kj(context, {prepared.force_group_electrostatics}))
    vdw = kj_per_mol_to_kcal_per_mol(_energy_kj(context, {prepared.force_group_vdw}))
    return elec, vdw


def _target_energy_delta_kcal(context, config: FEPConfig, lambda1: float, lambda2: float) -> float:
    set_context_lambda(context, config, lambda1)
    e1 = _energy_kj(context)
    set_context_lambda(context, config, lambda2)
    e2 = _energy_kj(context)
    return kj_per_mol_to_kcal_per_mol(e2 - e1)


def _window_energy_components(
    context, prepared: PreparedFEPSystem, lambda1: float, lambda2: float
) -> tuple[float, float, float, float, float]:
    config = prepared.config
    set_context_lambda(context, config, lambda1)
    e1 = _energy_kj(context)
    elec1, vdw1 = _component_energies_kcal(context, prepared)
    set_context_lambda(context, config, lambda2)
    e2 = _energy_kj(context)
    elec2, vdw2 = _component_energies_kcal(context, prepared)
    de = kj_per_mol_to_kcal_per_mol(e2 - e1)
    return elec1, elec2, vdw1, vdw2, de


def _instant_temperature_kelvin(simulation, prepared: PreparedFEPSystem) -> float:
    openmm = _require_openmm()
    system = prepared.system
    state = simulation.context.getState(getEnergy=True)
    kinetic = state.getKineticEnergy().value_in_unit(openmm.unit.kilojoule_per_mole)

    dof = 0
    for particle_index in range(system.getNumParticles()):
        mass = system.getParticleMass(particle_index).value_in_unit(openmm.unit.dalton)
        if mass > 0.0:
            dof += 3
    dof -= system.getNumConstraints()
    for force_index in range(system.getNumForces()):
        if isinstance(system.getForce(force_index), openmm.CMMotionRemover):
            dof -= 3

    if dof <= 0:
        return prepared.config.temperature
    gas_constant = openmm.unit.MOLAR_GAS_CONSTANT_R.value_in_unit(
        openmm.unit.kilojoule_per_mole / openmm.unit.kelvin
    )
    return 2.0 * kinetic / (dof * gas_constant)


SimulationFactory = Callable[[object], object]


def run_fep(prepared: PreparedFEPSystem, simulation_factory: SimulationFactory) -> None:
    """Run windows sequentially and write NAMD-like ``.fepout`` data.

    ``simulation_factory`` receives the prepared OpenMM System and must return an object with a
    ``context`` attribute and ``step(n)`` method, matching ``openmm.app.Simulation``.
    """

    simulation = simulation_factory(prepared.system)
    config = prepared.config
    with FEPOutWriter(config.alch_out_file, config.temperature) as writer:
        completed_steps = 0
        for window in config.iter_windows():
            _run_window(simulation, prepared, window, writer, completed_steps)
            completed_steps += window.steps


def _run_window(
    simulation,
    prepared: PreparedFEPSystem,
    window: LambdaWindow,
    writer: FEPOutWriter,
    window_start_step: int,
) -> None:
    config = prepared.config
    stats = FEPStatistics(config.temperature)
    back_stats = FEPStatistics(config.temperature)
    writer.start_window(window.lambda1, window.lambda2, window.lambda_idws)
    set_context_lambda(simulation.context, config, window.lambda1)

    local_step = 0
    if config.alch_equil_steps:
        equil_steps = min(config.alch_equil_steps, window.steps)
        if equil_steps:
            simulation.step(equil_steps)
            local_step = equil_steps
            stats = FEPStatistics(config.temperature)
            back_stats = FEPStatistics(config.temperature)
            if local_step == config.alch_equil_steps:
                writer.equilibration_done(config.alch_equil_steps, window.lambda1)
            set_context_lambda(simulation.context, config, window.lambda1)

    while local_step < window.steps:
        next_output_step = ((local_step // config.alch_out_freq) + 1) * config.alch_out_freq
        target_step = min(next_output_step, window.steps)
        simulation.step(target_step - local_step)
        local_step = target_step
        if local_step % config.alch_out_freq != 0:
            set_context_lambda(simulation.context, config, window.lambda1)
            continue

        absolute_step = window_start_step + local_step
        elec1, elec2, vdw1, vdw2, de = _window_energy_components(
            simulation.context, prepared, window.lambda1, window.lambda2
        )
        stats.add(de)
        set_context_lambda(simulation.context, config, window.lambda1)
        display_temperature = _instant_temperature_kelvin(simulation, prepared)
        writer.write_energy(
            absolute_step, elec1, elec2, vdw1, vdw2, de, stats, temperature=display_temperature
        )
        if window.lambda_idws is not None and window.sample_backward:
            elec1, elecb, vdw1, vdwb, de_back = _window_energy_components(
                simulation.context, prepared, window.lambda1, window.lambda_idws
            )
            back_stats.add(de_back)
            set_context_lambda(simulation.context, config, window.lambda1)
            writer.write_energy(
                absolute_step,
                elec1,
                elecb,
                vdw1,
                vdwb,
                de_back,
                back_stats,
                backward=True,
                temperature=display_temperature,
            )
        set_context_lambda(simulation.context, config, window.lambda1)
    writer.finish_window(window.lambda1, window.lambda2, stats.dG)
