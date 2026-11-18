// Adapted from Colvars PR 979. See ../UPSTREAM.md and ../LICENSE.
#include "colvarproxy_openmm.h"
#include "colvaratoms.h"
#include "colvar.h"
#include "colvarbias.h"
#include "openmm/Context.h"
#include "openmm/Integrator.h"
#include "openmm/OpenMMException.h"
#include "openmm/System.h"
#include "openmm/internal/ContextImpl.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
using namespace OpenMM;
using namespace ColvarsPlugin;
namespace {
std::recursive_mutex& colvarsMutex() {
    // Keep the mutex valid during Python/shared-library teardown.
    static auto* mutex = new std::recursive_mutex;
    return *mutex;
}
bool starts_with_record(const std::string& line)
{
    return line.rfind("ATOM  ", 0) == 0 || line.rfind("HETATM", 0) == 0;
}

double parse_pdb_double(const std::string& line, size_t offset, size_t count)
{
    if (line.size() < offset + count) {
        throw OpenMMException("PDB record is shorter than expected");
    }
    return std::stod(line.substr(offset, count));
}

char normalized_pdb_field(const std::string& field)
{
    if (field.empty()) {
        return '\0';
    }
    return static_cast<char>(std::toupper(static_cast<unsigned char>(field[0])));
}

double pdb_field_value(const std::string& line, char field)
{
    switch (field) {
        case 'O':
            return parse_pdb_double(line, 54, 6);
        case 'B':
            return parse_pdb_double(line, 60, 6);
        case 'X':
            return parse_pdb_double(line, 30, 8);
        case 'Y':
            return parse_pdb_double(line, 38, 8);
        case 'Z':
            return parse_pdb_double(line, 46, 8);
        default:
            throw OpenMMException("Unsupported PDB selector field");
    }
}

bool pdb_field_matches(const std::string& line, char field, double requestedValue)
{
    if (field == '\0') {
        return true;
    }
    const double value = pdb_field_value(line, field);
    if (requestedValue != 0.0) {
        return value == requestedValue;
    }
    return value != 0.0;
}

cvm::atom_pos pdb_position_angstrom(const std::string& line)
{
    return cvm::atom_pos(parse_pdb_double(line, 30, 8),
                         parse_pdb_double(line, 38, 8),
                         parse_pdb_double(line, 46, 8));
}

}

OpenMMColvarsLock::OpenMMColvarsLock(colvarproxy* active) : lock_(colvarsMutex()) {
    if (active) colvarmodule::proxy_static = active;
}

colvarproxy_openmm::colvarproxy_openmm(const ColvarsForce& owner, const System& system)
    : owner_(owner), system_(system),
      rng_(owner.getRandomSeed() == 0 ? std::random_device{}() : owner.getRandomSeed()) {
    OpenMMColvarsLock lock(this);
    version_int = get_version_from_string("2026-11-18");
    engine_name_ = "OpenMM";
    cvmodule = new colvarmodule(this);
    set_output_prefix(owner.getOutputPrefix());
    set_restart_output_prefix(owner.getOutputPrefix());
    if (!owner.getRestartInputPrefix().empty()) set_input_prefix(owner.getRestartInputPrefix());
    set_default_restart_frequency(0);
    add_config(owner.getConfigSource() == ColvarsForce::ConfigSource::File ? "configfile" : "config", owner.getConfig());
    usesPbc_ = owner.getPbcMode() == ColvarsForce::PbcMode::Enabled;
    if (owner.getPbcMode() == ColvarsForce::PbcMode::Auto) {
        for (int i = 0; i < system.getNumForces(); ++i) {
            const Force& f = system.getForce(i);
            if (&f != &owner && f.usesPeriodicBoundaryConditions()) usesPbc_ = true;
        }
    }
}

colvarproxy_openmm::~colvarproxy_openmm() {
    OpenMMColvarsLock lock(this);
    // Delete under the same lock used for evaluation. The base proxy destructor
    // skips this now-null pointer, preventing concurrent feature-table teardown.
    delete cvmodule;
    cvmodule = nullptr;
}

void colvarproxy_openmm::initializeColvars() {
    OpenMMColvarsLock lock(this);
    if (initialized_) return;
    if (owner_.getNumParticleCharges() != 0 && owner_.getNumParticleCharges() != system_.getNumParticles())
        throw OpenMMException("Colvars particle charges must contain exactly one entry per System particle");
    int result = set_unit_system(owner_.getUnitSystem(), false);
    result |= set_target_temperature(owner_.getTargetTemperature());
    result |= parse_module_config();
    result |= setup();
    result |= cvmodule->setup_input();
    result |= cvmodule->setup_output();
    if (result || cvmodule->get_error()) throw OpenMMException("Colvars initialization: " + get_error_msgs());
    // State mutation and fresh force evaluation must be separable for OpenMM.
    for (auto* cv : *cvmodule->variables()) {
        if (cv->get_time_step_factor() != 1)
            throw OpenMMException("OpenMM Colvars requires timeStepFactor 1 on every colvar; use OpenMM integrator force groups for multiple time stepping");
        if (cv->is_enabled(colvardeps::f_cv_extended_Lagrangian))
            throw OpenMMException("OpenMM Colvars: extendedLagrangian is not supported by repeat-safe evaluation");
    }
    if (cvmodule->scripted_forces())
        throw OpenMMException("OpenMM Colvars does not support scripted forces");
    for (auto* bias : cvmodule->biases) {
        if (bias->get_time_step_factor() != 1)
            throw OpenMMException("OpenMM Colvars requires timeStepFactor 1 on every bias; use OpenMM integrator force groups for multiple time stepping");
        if (!bias->supports_state_free_evaluation())
            throw OpenMMException("OpenMM Colvars bias does not support repeat-safe evaluation: " + bias->bias_type);
    }
    restartPending_ = cvmodule->has_loaded_state();
    initialized_ = true;
}

int colvarproxy_openmm::setup() {
    OpenMMColvarsLock lock(this);
    int result = colvarproxy::setup();
    return result | cvmodule->update_engine_parameters();
}

int colvarproxy_openmm::init_atom(int atomNumber) {
    const int id = check_atom_id(atomNumber);
    if (id < 0) return -1;
    auto found = std::find(atoms_ids.begin(), atoms_ids.end(), id);
    if (found != atoms_ids.end()) {
        const int index = static_cast<int>(found - atoms_ids.begin());
        ++atoms_refcount[index];
        return index;
    }
    const int index = add_atom_slot(id);
    atoms_masses[index] = system_.getParticleMass(id);
    atoms_charges[index] = owner_.getParticleCharge(id);
    return index;
}

int colvarproxy_openmm::check_atom_id(int atomNumber) {
    const int id = atomNumber - 1;
    if (id < 0 || id >= system_.getNumParticles())
        throw OpenMMException("Colvars atomNumbers are one-based: invalid atom " + std::to_string(atomNumber));
    return id;
}

int colvarproxy_openmm::request_atom_charges() {
    if (owner_.getNumParticleCharges() != system_.getNumParticles())
        throw OpenMMException("Charge-dependent Colvars components require setParticleCharges() before Context creation, with one fixed charge in elementary-charge units per System particle");
    return COLVARS_OK;
}

int colvarproxy_openmm::set_unit_system(const std::string& name, bool checkOnly) {
    if (name != "gromacs") {
        if (checkOnly) return COLVARS_ERROR;
        throw OpenMMException("OpenMM Colvars requires units gromacs (nm, kJ/mol); convert real-unit configurations explicitly");
    }
    if (!checkOnly) {
        units = "gromacs";
        angstrom_value_ = 0.1;
        kcal_mol_value_ = 4.184;
        boltzmann_ = 0.00831446261815324;
    }
    return COLVARS_OK;
}

void colvarproxy_openmm::request_total_force(bool yesno) {
    if (yesno) throw OpenMMException("OpenMM Colvars: total-force sampling (including ABF/eABF) is not supported; completed physical forces are not available in this evaluation path");
}

void colvarproxy_openmm::log(const std::string& message) {
    std::istringstream lines(message);
    std::string line;
    while (std::getline(lines, line)) std::cout << "colvars: " << line << '\n';
}
void colvarproxy_openmm::error(const std::string& message) {
    add_error_msg(message);
    throw OpenMMException("Colvars-OpenMM: " + message);
}
cvm::real colvarproxy_openmm::rand_gaussian() { return gaussian_(rng_); }

void colvarproxy_openmm::beginEvaluation(ContextImpl& context) {
    OpenMMColvarsLock lock(this);
    if (!initialized_) throw OpenMMException("Colvars proxy was not initialized");
    currentContext_ = &context;
    currentStep_ = context.getOwner().getStepCount();
    if (haveStep_ && currentStep_ < lastUpdatedStep_)
        throw OpenMMException("OpenMM Colvars cannot move backwards in step count; recreate Context and load a Colvars restart");
    if (!haveStep_) {
        initialOpenmmStep_ = currentStep_;
        initialColvarsStep_ = restartPending_ ? cvmodule->it_restart : currentStep_;
    }
    cvmodule->it = initialColvarsStep_ + (currentStep_ - initialOpenmmStep_);
    const double timestepFs = 1000.0 * context.getIntegrator().getStepSize();
    if (timestepFs != lastTimestepFs_) {
        set_integration_timestep(timestepFs);
        if (cvmodule->update_engine_parameters())
            throw OpenMMException("Colvars could not update the integration timestep");
        lastTimestepFs_ = timestepFs;
    }
    if (usesPbc_) {
        Vec3 a,b,c;
        context.getPeriodicBoxVectors(a,b,c);
        boundaries_.set_boundaries(1,1,1, cvm::rvector(a[0],a[1],a[2]),
                                    cvm::rvector(b[0],b[1],b[2]),
                                    cvm::rvector(c[0],c[1],c[2]));
    } else boundaries_.reset();
    biasEnergy_ = 0.0;
    auto& forces = *modify_atom_applied_forces();
    std::fill(forces.begin(), forces.end(), cvm::rvector(0.0,0.0,0.0));
}

double colvarproxy_openmm::evaluateColvars() {
    OpenMMColvarsLock lock(this);
    // The step saved in a restart has already updated its adaptive history.
    // Its first evaluation restores current forces without depositing again.
    const bool advance = haveStep_ ? currentStep_ != lastUpdatedStep_ : !restartPending_;
    const int result = cvmodule->calc(advance);
    if (result || cvmodule->get_error()) throw OpenMMException("Colvars evaluation: " + get_error_msgs());
    haveStep_ = true;
    restartPending_ = false;
    lastUpdatedStep_ = currentStep_;
    if (advance) ++biasUpdateCount_;
    return biasEnergy_;
}

double colvarproxy_openmm::calculate(ContextImpl& context, const std::vector<Vec3>& positions) {
    OpenMMColvarsLock lock(this);
    beginEvaluation(context);
    auto& selected = *modify_atom_positions();
    for (size_t i = 0; i < atoms_ids.size(); ++i) {
        const Vec3& p = positions.at(atoms_ids[i]);
        selected[i] = cvm::rvector(p[0],p[1],p[2]);
    }
    return evaluateColvars();
}

void colvarproxy_openmm::addForcesToOpenMM(std::vector<Vec3>& forces) const {
    const auto& applied = *get_atom_applied_forces();
    for (size_t i = 0; i < atoms_ids.size(); ++i)
        forces.at(atoms_ids[i]) += Vec3(applied[i].x,applied[i].y,applied[i].z);
}

double colvarproxy_openmm::getColvarValue(const std::string& name) const {
    OpenMMColvarsLock lock(const_cast<colvarproxy_openmm*>(this));
    if (!haveStep_) throw OpenMMException("Evaluate energy or forces before reading a Colvars value");
    auto* cv = cvmodule->colvar_by_name(name);
    if (!cv) throw OpenMMException("Unknown Colvars variable: " + name);
    if (cv->value().type() != colvarvalue::type_scalar)
        throw OpenMMException("getColvarValue requires a scalar variable");
    return cv->value().real_value;
}
std::string colvarproxy_openmm::getBiasState() {
    OpenMMColvarsLock lock(this);
    std::string result;
    if (cvmodule->write_restart_string(result)) throw OpenMMException("Cannot serialize Colvars state");
    return result;
}
void colvarproxy_openmm::flushOutputs() {
    OpenMMColvarsLock lock(this);
    if (!haveStep_) throw OpenMMException("Evaluate energy or forces before flushing Colvars output");
    // Use the core end-of-run output path explicitly, while exceptions can be
    // reported to Python. Destructors must not perform fallible file I/O.
    if (post_run() || cvmodule->get_error())
        throw OpenMMException("Cannot flush Colvars output: " + get_error_msgs());
}
long long colvarproxy_openmm::getBiasUpdateCount() const {
    OpenMMColvarsLock lock(const_cast<colvarproxy_openmm*>(this));
    return biasUpdateCount_;
}

int colvarproxy_openmm::load_atoms_pdb(char const* filename,
                                      cvm::atom_group& atoms,
                                      std::string const& pdbField,
                                      double pdbFieldValue)
{
    const char field = normalized_pdb_field(pdbField);
    if (field == '\0') {
        return cvmodule->error("Error: load_atoms_pdb requires a PDB selector field.\n",
                               COLVARS_INPUT_ERROR);
    }

    const std::string filenameString = filename == nullptr ? std::string() : std::string(filename);
    std::ifstream input(filenameString);
    if (!input) {
        return cvmodule->error("Error: cannot open PDB file \"" + filenameString + "\".\n",
                               COLVARS_FILE_ERROR);
    }

    auto modifyAtoms = atoms.get_atom_modifier();
    std::string line;
    int pdbAtomIndex = 0;
    while (std::getline(input, line)) {
        if (!starts_with_record(line)) {
            continue;
        }
        if (pdb_field_matches(line, field, pdbFieldValue)) {
            if (pdbAtomIndex >= system_.getNumParticles()) {
                return cvmodule->error(
                    "Error: PDB file \"" + filenameString + "\" selects atom number " +
                        std::to_string(pdbAtomIndex + 1) +
                        ", which is outside the OpenMM System.\n",
                    COLVARS_INPUT_ERROR);
            }
            if (atoms.is_enabled(colvardeps::f_ag_scalable)) {
                modifyAtoms.add_atom_id(pdbAtomIndex);
            } else {
                modifyAtoms.add_atom(cvm::atom_group::init_atom_from_proxy(this, pdbAtomIndex + 1));
            }
        }
        ++pdbAtomIndex;
    }
    return cvmodule->get_error() ? COLVARS_ERROR : COLVARS_OK;
}

int colvarproxy_openmm::load_coords_pdb(char const* filename,
                                       std::vector<cvm::atom_pos>& positions,
                                       std::vector<int> const& sortedIds,
                                       std::string const& pdbField,
                                       double pdbFieldValue)
{
    const char field = normalized_pdb_field(pdbField);
    const std::string filenameString = filename == nullptr ? std::string() : std::string(filename);
    int previousId = -1;
    for (int atomId : sortedIds) {
        if (atomId < 0 || atomId >= system_.getNumParticles()) {
            return cvmodule->error(
                "Error: PDB file \"" + filenameString + "\" requested atom number " +
                    std::to_string(atomId + 1) + ", which is outside the OpenMM System.\n",
                COLVARS_INPUT_ERROR);
        }
        if (atomId <= previousId) {
            return cvmodule->error(
                "Error: requested atom numbers for PDB file \"" + filenameString +
                    "\" must be strictly increasing without duplicates.\n",
                COLVARS_INPUT_ERROR);
        }
        previousId = atomId;
    }

    std::ifstream input(filenameString);
    if (!input) {
        return cvmodule->error("Error: cannot open PDB file \"" + filenameString + "\".\n",
                               COLVARS_FILE_ERROR);
    }

    std::vector<cvm::atom_pos> selectedPositions;
    std::string line;
    int pdbAtomIndex = 0;
    size_t nextSortedId = 0;
    while (std::getline(input, line)) {
        if (!starts_with_record(line)) {
            continue;
        }

        bool useThisAtom = pdb_field_matches(line, field, pdbFieldValue);
        if (!sortedIds.empty()) {
            useThisAtom = nextSortedId < sortedIds.size() && sortedIds[nextSortedId] == pdbAtomIndex;
            if (useThisAtom) {
                ++nextSortedId;
            }
        }

        if (useThisAtom) {
            selectedPositions.push_back(angstrom_value_ * pdb_position_angstrom(line));
        }
        ++pdbAtomIndex;
    }

    if (!sortedIds.empty() && nextSortedId != sortedIds.size()) {
        return cvmodule->error("Error: PDB file \"" + filenameString +
                                   "\" does not contain requested atom number " +
                                   std::to_string(sortedIds[nextSortedId] + 1) + ".\n",
                               COLVARS_INPUT_ERROR);
    }
    if (!positions.empty() && positions.size() != selectedPositions.size()) {
        return cvmodule->error("Error: PDB coordinate count does not match requested positions.\n",
                               COLVARS_INPUT_ERROR);
    }

    positions = selectedPositions;
    return COLVARS_OK;
}
