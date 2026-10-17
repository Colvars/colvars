#include "openmm/internal/OpenMMColvarProxy.h"

#include "colvars_units.h"

#include "openmm/Context.h"
#include "openmm/CustomIntegrator.h"
#include "openmm/Force.h"
#include "openmm/Integrator.h"
#include "openmm/OpenMMException.h"
#include "openmm/System.h"
#include "openmm/internal/ContextImpl.h"

#include "colvaratoms.h"
#include "colvarmodule.h"
#include "colvars_version.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>

using namespace OpenMM;

namespace {

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

bool is_orthorhombic(const Vec3& a, const Vec3& b, const Vec3& c)
{
    constexpr double kTolerance = 1.0e-12;
    return std::abs(a[1]) < kTolerance && std::abs(a[2]) < kTolerance &&
           std::abs(b[0]) < kTolerance && std::abs(b[2]) < kTolerance &&
           std::abs(c[0]) < kTolerance && std::abs(c[1]) < kTolerance;
}

bool same_vec3(const Vec3& lhs, const Vec3& rhs)
{
    constexpr double kTolerance = 1.0e-12;
    return std::abs(lhs[0] - rhs[0]) < kTolerance &&
           std::abs(lhs[1] - rhs[1]) < kTolerance &&
           std::abs(lhs[2] - rhs[2]) < kTolerance;
}

bool same_position(const Vec3& lhs, const Vec3& rhs)
{
    // CUDA mixed precision can return the same stored position with sub-micro-nanometer
    // differences when it is downloaded again for a reporter force evaluation.
    constexpr double kToleranceNm = 1.0e-6;
    return std::abs(lhs[0] - rhs[0]) < kToleranceNm &&
           std::abs(lhs[1] - rhs[1]) < kToleranceNm &&
           std::abs(lhs[2] - rhs[2]) < kToleranceNm;
}

bool same_scalar(double lhs, double rhs)
{
    if (lhs == rhs) {
        return true;
    }
    if (!std::isfinite(lhs) || !std::isfinite(rhs)) {
        return false;
    }
    constexpr double kTolerance = 1.0e-12;
    const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
    return std::abs(lhs - rhs) < kTolerance * scale;
}

std::string format_colvars_error(const std::string& message)
{
    if (message.empty()) {
        return "Error in Colvars-OpenMM interface";
    }
    return "Colvars-OpenMM: " + message;
}

bool is_repeated_extended_coordinate_message(const std::string& line)
{
    return line.find("Reinitializing extended coordinate to last value") != std::string::npos ||
           line.find("Reinitializing extended coordinate to colvar value") != std::string::npos;
}

std::uint64_t splitmix64(std::uint64_t& state)
{
    state += 0x9e3779b97f4a7c15ULL;
    std::uint64_t value = state;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

std::uint64_t make_random_seed(int seed)
{
    if (seed != 0) {
        return static_cast<std::uint64_t>(seed);
    }
    const auto now = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::uint64_t state = static_cast<std::uint64_t>(now);
    return splitmix64(state);
}

cvm::real uniform_open_0_1(std::uint64_t& state)
{
    const std::uint64_t value = splitmix64(state);
    const double uniform = static_cast<double>(value >> 11) * (1.0 / 9007199254740992.0);
    return static_cast<cvm::real>(std::max(uniform, std::numeric_limits<double>::min()));
}

std::string strip_colvars_comments(const std::string& input)
{
    std::string output;
    output.reserve(input.size());
    bool inComment = false;
    for (char ch : input) {
        if (inComment) {
            if (ch == '\n' || ch == '\r') {
                inComment = false;
                output.push_back(' ');
            }
            continue;
        }
        if (ch == '#') {
            inComment = true;
            continue;
        }
        output.push_back(ch);
    }
    return output;
}

std::vector<std::string> tokenize_colvars_config(const std::string& input)
{
    std::vector<std::string> tokens;
    std::string token;
    for (unsigned char ch : input) {
        if (ch == '{' || ch == '}') {
            if (!token.empty()) {
                tokens.push_back(token);
                token.clear();
            }
            tokens.emplace_back(1, static_cast<char>(ch));
            continue;
        }
        if (std::isspace(ch)) {
            if (!token.empty()) {
                tokens.push_back(token);
                token.clear();
            }
            continue;
        }
        token.push_back(static_cast<char>(std::tolower(ch)));
    }
    if (!token.empty()) {
        tokens.push_back(token);
    }
    return tokens;
}

std::vector<std::vector<std::string>> top_level_blocks(const std::vector<std::string>& tokens,
                                                       const std::string& keyword)
{
    std::vector<std::vector<std::string>> blocks;
    for (size_t index = 0; index + 1 < tokens.size(); ++index) {
        if (tokens[index] != keyword || tokens[index + 1] != "{") {
            continue;
        }
        int depth = 1;
        const size_t begin = index + 2;
        size_t end = begin;
        for (; end < tokens.size(); ++end) {
            if (tokens[end] == "{") {
                ++depth;
            } else if (tokens[end] == "}") {
                --depth;
                if (depth == 0) {
                    blocks.emplace_back(tokens.begin() + static_cast<std::ptrdiff_t>(begin),
                                        tokens.begin() + static_cast<std::ptrdiff_t>(end));
                    index = end;
                    break;
                }
            }
        }
    }
    return blocks;
}

bool token_value_is_on(const std::vector<std::string>& tokens, const std::string& key)
{
    for (size_t index = 0; index + 1 < tokens.size(); ++index) {
        if (tokens[index] == key) {
            return tokens[index + 1] == "on" || tokens[index + 1] == "yes" ||
                   tokens[index + 1] == "true" || tokens[index + 1] == "1";
        }
    }
    return false;
}

std::string token_value_after(const std::vector<std::string>& tokens, const std::string& key)
{
    for (size_t index = 0; index + 1 < tokens.size(); ++index) {
        if (tokens[index] == key && tokens[index + 1] != "{") {
            return tokens[index + 1];
        }
    }
    return {};
}

bool contains_token(const std::vector<std::string>& tokens, const std::string& value)
{
    return std::find(tokens.begin(), tokens.end(), value) != tokens.end();
}

std::vector<std::string> abf_colvars(const std::vector<std::string>& tokens)
{
    static const std::set<std::string> abfKeys = {
        "name", "stepzerodata", "timestepfactor", "outputenergy", "outputfreq",
        "applybias", "hidejacobian", "fullsamples", "minsamples", "inputprefix",
        "historyfreq", "shared", "sharedfreq", "updatebias", "maxforce",
        "czarestimator", "integrate", "uiestimator"
    };

    for (size_t index = 0; index + 1 < tokens.size(); ++index) {
        if (tokens[index] != "colvars") {
            continue;
        }
        std::vector<std::string> names;
        size_t cursor = index + 1;
        if (tokens[cursor] == "{") {
            ++cursor;
            while (cursor < tokens.size() && tokens[cursor] != "}") {
                if (tokens[cursor] != "{") {
                    names.push_back(tokens[cursor]);
                }
                ++cursor;
            }
            return names;
        }
        while (cursor < tokens.size() && tokens[cursor] != "}" &&
               abfKeys.find(tokens[cursor]) == abfKeys.end()) {
            if (tokens[cursor] != "{") {
                names.push_back(tokens[cursor]);
            }
            ++cursor;
        }
        return names;
    }
    return {};
}

bool config_has_geometrical_eabf(const std::string& configText)
{
    const auto tokens = tokenize_colvars_config(strip_colvars_comments(configText));
    std::set<std::string> geometricalExtendedColvars;
    int unnamedColvarIndex = 0;

    for (const auto& block : top_level_blocks(tokens, "colvar")) {
        std::string name = token_value_after(block, "name");
        if (name.empty()) {
            name = "colvar" + std::to_string(unnamedColvarIndex);
        }
        ++unnamedColvarIndex;

        const bool extended = token_value_is_on(block, "extendedlagrangian");
        const bool external = contains_token(block, "alchlambda") || contains_token(block, "tabfalpha");
        if (extended && !external) {
            geometricalExtendedColvars.insert(name);
        }
    }

    if (geometricalExtendedColvars.empty()) {
        return false;
    }

    for (const auto& block : top_level_blocks(tokens, "abf")) {
        for (const std::string& name : abf_colvars(block)) {
            if (geometricalExtendedColvars.find(name) != geometricalExtendedColvars.end()) {
                return true;
            }
        }
    }
    return false;
}

std::string owner_config_text(const ColvarsForce& owner)
{
    if (owner.getConfigSource() == ColvarsForce::ConfigSource::String) {
        return owner.getConfig();
    }
    std::ifstream configFile(owner.getConfig());
    if (!configFile) {
        return {};
    }
    std::ostringstream buffer;
    buffer << configFile.rdbuf();
    return buffer.str();
}

bool config_tabf_alpha_boundaries(const std::string& configText,
                                  cvm::real& lower,
                                  cvm::real& upper)
{
    const auto tokens = tokenize_colvars_config(strip_colvars_comments(configText));
    for (const auto& block : top_level_blocks(tokens, "colvar")) {
        if (!contains_token(block, "tabfalpha")) {
            continue;
        }
        const std::string lowerText = token_value_after(block, "lowerboundary");
        const std::string upperText = token_value_after(block, "upperboundary");
        if (lowerText.empty() || upperText.empty()) {
            return false;
        }
        try {
            const double parsedLower = std::stod(lowerText);
            const double parsedUpper = std::stod(upperText);
            if (!std::isfinite(parsedLower) || !std::isfinite(parsedUpper) ||
                parsedLower >= parsedUpper) {
                return false;
            }
            lower = static_cast<cvm::real>(parsedLower);
            upper = static_cast<cvm::real>(parsedUpper);
            return true;
        } catch (const std::exception&) {
            return false;
        }
    }
    return false;
}

} // namespace

OpenMMColvarProxy::OpenMMColvarProxy(const ColvarsForce& owner, const System& system)
    : owner_(owner),
      system_(system),
      currentContext_(nullptr),
      biasEnergyKcalPerMol_(0.0),
      isInitialized_(false),
      autoPbcResolved_(false),
      cachedAutoPbc_(false),
      pbcLatticeInitialized_(false),
      reportedRepeatedExtendedCoordinate_(false),
      totalForcesSameStep_(true),
      stepMappingInitialized_(false),
      openmmInitialStep_(0),
      colvarsInitialStep_(0),
      lastLambdaWritebackStep_(-1),
      lastTabfAlphaWritebackStep_(-1),
      lastIntegrationTimestepFs_(-1.0),
      gaussianRngState_(make_random_seed(owner_.getRandomSeed())),
      gaussianRngHasSpare_(false),
      gaussianRngSpare_(0.0),
      currentLambdaDerivativeValid_(false),
      currentLambdaDerivativeKJPerMol_(0.0),
      tabfAlphaBoundariesValid_(false),
      tabfAlphaLowerBoundary_(0.0),
      tabfAlphaUpperBoundary_(0.0),
      cacheValid_(false),
      cachedOpenmmStep_(-1),
      cachedUsesPbc_(false),
      cachedLambdaValue_(0.0),
      cachedTabfAlphaValue_(1.0),
      cachedBiasEnergyKJPerMol_(0.0)
{
    tabfAlphaBoundariesValid_ = config_tabf_alpha_boundaries(
        owner_config_text(owner_), tabfAlphaLowerBoundary_, tabfAlphaUpperBoundary_);

    version_int = get_version_from_string(COLVARS_VERSION);
    engine_name_ = "OpenMM";
    engine_ready_ = true;
    b_simulation_running = true;
    b_simulation_continuing = false;
    use_internal_pbc_ = true;
    updated_masses_ = true;
    updated_charges_ = true;

    cvmodule = new colvarmodule(this);
    cvmodule->log("Using OpenMM interface.\n");
    totalForcesSameStep_ = inferTotalForcesSameStep();

    colvarproxy_io::set_output_prefix(owner_.getOutputPrefix());
    colvarproxy_io::set_restart_output_prefix(owner_.getRestartOutputPrefix());
    if (!owner_.getRestartInputPrefix().empty()) {
        colvarproxy_io::set_input_prefix(owner_.getRestartInputPrefix());
    }
    colvarproxy_io::set_default_restart_frequency(owner_.getRestartFrequency());

    if (owner_.getConfigSource() == ColvarsForce::ConfigSource::File) {
        add_config("configfile", owner_.getConfig());
    } else {
        add_config("config", owner_.getConfig());
    }

    initializeColvars();
}

OpenMMColvarProxy::~OpenMMColvarProxy()
{
    close_output_streams();
}

void OpenMMColvarProxy::initializeColvars()
{
    int errorCode = COLVARS_OK;
    errorCode |= set_unit_system(owner_.getUnitSystem(), false);
    errorCode |= set_target_temperature(owner_.getTargetTemperature());
    errorCode |= parse_module_config();
    errorCode |= setup();
    errorCode |= cvmodule->setup_input();
    errorCode |= cvmodule->setup_output();
    if (errorCode != COLVARS_OK || cvmodule->get_error() != COLVARS_OK) {
        throw OpenMMException(format_colvars_error(get_error_msgs()));
    }
    isInitialized_ = true;
}

int OpenMMColvarProxy::setup()
{
    int errorCode = colvarproxy::setup();
    errorCode |= cvmodule->update_engine_parameters();
    return errorCode;
}

int OpenMMColvarProxy::init_atom(int atomNumber)
{
    int atomId = atomNumber - 1;
    const auto cachedIndex = atomIndexById_.find(atomId);
    if (cachedIndex != atomIndexById_.end()) {
        atoms_refcount[cachedIndex->second] += 1;
        return cachedIndex->second;
    }

    atomId = check_atom_id(atomNumber);
    if (atomId < 0 || cvmodule->get_error() != COLVARS_OK) {
        return COLVARS_INPUT_ERROR;
    }

    const int proxyIndex = add_atom_slot(atomId);
    atomIndexById_[atomId] = proxyIndex;
    updateAtomProperties(proxyIndex);
    invalidateCache();
    return proxyIndex;
}

void OpenMMColvarProxy::clear_atom(int index)
{
    colvarproxy::clear_atom(index);
    updateAtomIndexMap();
    invalidateCache();
}

int OpenMMColvarProxy::check_atom_id(int atomNumber)
{
    const int atomId = atomNumber - 1;
    if (atomId < 0 || atomId >= system_.getNumParticles()) {
        cvmodule->error("Error: invalid atom number specified, " + cvm::to_str(atomNumber) + "\n",
                        COLVARS_INPUT_ERROR);
        return -1;
    }
    return atomId;
}

void OpenMMColvarProxy::request_total_force(bool yesno)
{
    total_force_requested = yesno;
}

bool OpenMMColvarProxy::total_forces_enabled() const
{
    return total_force_requested;
}

bool OpenMMColvarProxy::total_forces_same_step() const
{
    return totalForcesSameStep_;
}

int OpenMMColvarProxy::set_unit_system(std::string const& unitsIn, bool checkOnly)
{
    if (checkOnly) {
        if ((units != "" && unitsIn != units) || (units == "" && unitsIn != "real")) {
            cvmodule->error("Specified unit system \"" + unitsIn +
                            "\" is incompatible with the OpenMM Colvars interface.\n",
                            COLVARS_INPUT_ERROR);
            return COLVARS_ERROR;
        }
        return COLVARS_OK;
    }

    if (unitsIn != "real") {
        cvmodule->error("Error: the OpenMM Colvars v1 interface currently supports only "
                        "\"real\" units (A, kcal/mol).\n",
                        COLVARS_INPUT_ERROR);
        return COLVARS_ERROR;
    }

    units = unitsIn;
    angstrom_value_ = 1.0;
    kcal_mol_value_ = 1.0;
    boltzmann_ = 0.001987191;
    return COLVARS_OK;
}

void OpenMMColvarProxy::add_energy(cvm::real energy)
{
    biasEnergyKcalPerMol_ += energy;
}

void OpenMMColvarProxy::log(std::string const& message)
{
    std::istringstream lines(message);
    std::string line;
    while (std::getline(lines, line)) {
        if (is_repeated_extended_coordinate_message(line)) {
            if (!reportedRepeatedExtendedCoordinate_) {
                std::cout << "colvars: Reinitializing extended coordinate after repeated "
                             "OpenMM force evaluation at the same step; further messages suppressed.\n";
                reportedRepeatedExtendedCoordinate_ = true;
            }
            continue;
        }
        std::cout << "colvars: " << line << '\n';
    }
}

void OpenMMColvarProxy::error(std::string const& message)
{
    add_error_msg(message);
    throw OpenMMException(format_colvars_error(message));
}

cvm::real OpenMMColvarProxy::rand_gaussian()
{
    if (gaussianRngHasSpare_) {
        gaussianRngHasSpare_ = false;
        return gaussianRngSpare_;
    }

    const cvm::real u1 = uniform_open_0_1(gaussianRngState_);
    const cvm::real u2 = uniform_open_0_1(gaussianRngState_);
    constexpr cvm::real twoPi = 6.283185307179586476925286766559;
    const cvm::real radius = std::sqrt(-2.0 * std::log(u1));
    const cvm::real angle = twoPi * u2;
    gaussianRngSpare_ = radius * std::sin(angle);
    gaussianRngHasSpare_ = true;
    return radius * std::cos(angle);
}

int OpenMMColvarProxy::get_alch_lambda(cvm::real* lambda)
{
    if (!owner_.getAlchemyEnabled()) {
        return cvmodule->error("Error: alchemical lambda support is disabled for this ColvarsForce.\n",
                               COLVARS_NOT_IMPLEMENTED);
    }
    if (lambda == nullptr) {
        return COLVARS_BUG_ERROR;
    }
    if (currentContext_ == nullptr) {
        *lambda = 0.0;
        return COLVARS_OK;
    }
    *lambda = contextParameter(*currentContext_, owner_.getLambdaParameterName(), 0.0);
    return COLVARS_OK;
}

int OpenMMColvarProxy::send_alch_lambda()
{
    if (!owner_.getAlchemyEnabled()) {
        return cvmodule->error("Error: alchemical lambda support is disabled for this ColvarsForce.\n",
                               COLVARS_NOT_IMPLEMENTED);
    }
    if (currentContext_ != nullptr && cached_alch_lambda_changed) {
        setContextParameter(*currentContext_, owner_.getLambdaParameterName(), cached_alch_lambda);
        cached_alch_lambda_changed = false;
    }
    return COLVARS_OK;
}

int OpenMMColvarProxy::request_alch_energy_freq(int freq)
{
    if (!owner_.getAlchemyEnabled()) {
        return cvmodule->error("Error: alchemical lambda support is disabled for this ColvarsForce.\n",
                               COLVARS_NOT_IMPLEMENTED);
    }
    if (freq != 1) {
        return cvmodule->error("Error: alchemical Colvars in OpenMM require timeStepFactor 1.\n",
                               COLVARS_NOT_IMPLEMENTED);
    }
    return COLVARS_OK;
}

int OpenMMColvarProxy::get_dE_dlambda(cvm::real* dEdLambda)
{
    if (!owner_.getAlchemyEnabled()) {
        return cvmodule->error("Error: alchemical lambda support is disabled for this ColvarsForce.\n",
                               COLVARS_NOT_IMPLEMENTED);
    }
    if (dEdLambda == nullptr) {
        return COLVARS_BUG_ERROR;
    }
    if (currentContext_ == nullptr) {
        *dEdLambda = 0.0;
        return COLVARS_OK;
    }

    double derivative_kj_per_mol = currentLambdaDerivativeKJPerMol_;
    if (!currentLambdaDerivativeValid_) {
        std::map<std::string, double> derivatives;
        currentContext_->getEnergyParameterDerivatives(derivatives);
        const auto derivative = derivatives.find(owner_.getLambdaParameterName());
        if (derivative != derivatives.end()) {
            derivative_kj_per_mol = derivative->second;
        } else {
            derivative_kj_per_mol =
                contextParameter(*currentContext_, owner_.getLambdaDerivativeParameterName(), 0.0);
        }
    }
    if (owner_.getScaleLambdaDerivativeByTabfAlpha()) {
        const double alpha = tabfAlpha(*currentContext_);
        if (alpha <= 0.0) {
            return cvmodule->error("Error: TABF alpha must be positive when scaling dE/dlambda.\n",
                                   COLVARS_ERROR);
        }
        derivative_kj_per_mol /= alpha;
    }
    *dEdLambda = colvars_units::kj_per_mol_to_kcal_per_mol(derivative_kj_per_mol);
    return COLVARS_OK;
}

#ifdef OPENMM_COLVARS_HAS_TABF_ALPHA
int OpenMMColvarProxy::get_tabf_alpha(cvm::real* alpha)
{
    if (!owner_.getTabfEnabled()) {
        return cvmodule->error("Error: TABF alpha support is disabled for this ColvarsForce.\n",
                               COLVARS_NOT_IMPLEMENTED);
    }
    if (alpha == nullptr) {
        return COLVARS_BUG_ERROR;
    }
    if (currentContext_ == nullptr) {
        *alpha = 1.0;
        return COLVARS_OK;
    }
    *alpha = tabfAlpha(*currentContext_);
    return COLVARS_OK;
}

#ifdef OPENMM_COLVARS_HAS_TABF_ALPHA_BOUNDARIES
int OpenMMColvarProxy::get_tabf_alpha_boundaries(cvm::real* lower, cvm::real* upper)
{
    if (!owner_.getTabfEnabled()) {
        return cvmodule->error("Error: TABF alpha support is disabled for this ColvarsForce.\n",
                               COLVARS_NOT_IMPLEMENTED);
    }
    if (lower == nullptr || upper == nullptr) {
        return COLVARS_BUG_ERROR;
    }
    if (!tabfAlphaBoundariesValid_) {
        return cvmodule->error(
            "Error: OpenMM TABF alpha requires finite lowerBoundary and upperBoundary "
            "values with lowerBoundary < upperBoundary in the tabfAlpha colvar.\n",
            COLVARS_INPUT_ERROR);
    }
    *lower = tabfAlphaLowerBoundary_;
    *upper = tabfAlphaUpperBoundary_;
    return COLVARS_OK;
}
#endif

int OpenMMColvarProxy::send_tabf_alpha()
{
    if (!owner_.getTabfEnabled()) {
        return cvmodule->error("Error: TABF alpha support is disabled for this ColvarsForce.\n",
                               COLVARS_NOT_IMPLEMENTED);
    }
    if (currentContext_ != nullptr && cached_tabf_alpha_changed) {
        if (cached_tabf_alpha <= 0.0) {
            return cvmodule->error("Error: TABF alpha must remain positive.\n", COLVARS_ERROR);
        }
        const bool wroteIntegrator =
            setIntegratorGlobal(*currentContext_, owner_.getTabfAlphaParameterName(), cached_tabf_alpha);
        if (!wroteIntegrator) {
            setContextParameter(*currentContext_, owner_.getTabfAlphaParameterName(), cached_tabf_alpha);
        } else {
            try {
                setContextParameter(*currentContext_, owner_.getTabfAlphaParameterName(), cached_tabf_alpha);
            } catch (const OpenMMException&) {
            }
        }
        cached_tabf_alpha_changed = false;
    }
    return COLVARS_OK;
}

int OpenMMColvarProxy::request_tabf_energy_freq(int freq)
{
    if (!owner_.getTabfEnabled()) {
        return cvmodule->error("Error: TABF alpha support is disabled for this ColvarsForce.\n",
                               COLVARS_NOT_IMPLEMENTED);
    }
    if (freq != 1) {
        return cvmodule->error("Error: TABF alpha variables in OpenMM require timeStepFactor 1.\n",
                               COLVARS_NOT_IMPLEMENTED);
    }
    return COLVARS_OK;
}

int OpenMMColvarProxy::get_dE_dtabf_alpha(cvm::real* dEdAlpha)
{
    if (!owner_.getTabfEnabled()) {
        return cvmodule->error("Error: TABF alpha support is disabled for this ColvarsForce.\n",
                               COLVARS_NOT_IMPLEMENTED);
    }
    if (dEdAlpha == nullptr) {
        return COLVARS_BUG_ERROR;
    }
    if (currentContext_ == nullptr) {
        *dEdAlpha = 0.0;
        return COLVARS_OK;
    }

    double derivative_kj_per_mol =
        integratorGlobal(*currentContext_, owner_.getTabfAlphaDerivativeParameterName(),
                         std::numeric_limits<double>::quiet_NaN());
    if (std::isnan(derivative_kj_per_mol)) {
        derivative_kj_per_mol =
            contextParameter(*currentContext_, owner_.getTabfAlphaDerivativeParameterName(), 0.0);
    }
    *dEdAlpha = colvars_units::kj_per_mol_to_kcal_per_mol(derivative_kj_per_mol);
    return COLVARS_OK;
}
#endif

int OpenMMColvarProxy::load_atoms_pdb(char const* filename,
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

int OpenMMColvarProxy::load_coords_pdb(char const* filename,
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
            selectedPositions.push_back(pdb_position_angstrom(line));
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

int OpenMMColvarProxy::check_replicas_enabled()
{
    return COLVARS_NOT_IMPLEMENTED;
}

int OpenMMColvarProxy::replica_index()
{
    return 0;
}

int OpenMMColvarProxy::num_replicas()
{
    return 1;
}

void OpenMMColvarProxy::replica_comm_barrier()
{
}

int OpenMMColvarProxy::replica_comm_recv(char*, int, int)
{
    return COLVARS_NOT_IMPLEMENTED;
}

int OpenMMColvarProxy::replica_comm_send(char*, int, int)
{
    return COLVARS_NOT_IMPLEMENTED;
}

int OpenMMColvarProxy::run_force_callback()
{
    return cvmodule->error("Error: scripted Colvars forces are not supported by the OpenMM v1 interface.\n",
                           COLVARS_NOT_IMPLEMENTED);
}

int OpenMMColvarProxy::run_colvar_callback(std::string const&,
                                           std::vector<const colvarvalue*> const&,
                                           colvarvalue&)
{
    return cvmodule->error("Error: scripted Colvars functions are not supported by the OpenMM v1 interface.\n",
                           COLVARS_NOT_IMPLEMENTED);
}

int OpenMMColvarProxy::run_colvar_gradient_callback(std::string const&,
                                                    std::vector<const colvarvalue*> const&,
                                                    std::vector<cvm::matrix2d<cvm::real>>&)
{
    return cvmodule->error("Error: scripted Colvars gradients are not supported by the OpenMM v1 interface.\n",
                           COLVARS_NOT_IMPLEMENTED);
}

int OpenMMColvarProxy::updateInput(ContextImpl& context, const std::vector<Vec3>& positions)
{
    currentContext_ = &context;
    biasEnergyKcalPerMol_ = 0.0;

    updateColvarsStep(context);
    const double timestepFs = 1000.0 * context.getIntegrator().getStepSize();
    if (timestepFs != lastIntegrationTimestepFs_) {
        set_integration_timestep(timestepFs);
        lastIntegrationTimestepFs_ = timestepFs;
    }
    updatePbc(context);
    updatePositions(positions);
    captureLambdaDerivative(context);
    updateTotalForces(context);

    auto& appliedForces = *(modify_atom_applied_forces());
    std::fill(appliedForces.begin(), appliedForces.end(), cvm::rvector(0.0, 0.0, 0.0));
    return COLVARS_OK;
}

double OpenMMColvarProxy::calculate(ContextImpl& context, const std::vector<Vec3>& positions)
{
    if (!isInitialized_) {
        initializeColvars();
    }

    if (currentContext_ != &context) {
        lastLambdaWritebackStep_ = -1;
        lastTabfAlphaWritebackStep_ = -1;
    }
    currentContext_ = &context;

    long long openmmStep = 0;
    bool usesPbc = false;
    Vec3 boxA_nm;
    Vec3 boxB_nm;
    Vec3 boxC_nm;
    double lambdaValue = 0.0;
    double tabfAlphaValue = 1.0;
    std::vector<Vec3> requestedPositions_nm;
    const bool canUseCache = buildCacheKey(context,
                                           positions,
                                           openmmStep,
                                           usesPbc,
                                           boxA_nm,
                                           boxB_nm,
                                           boxC_nm,
                                           lambdaValue,
                                           tabfAlphaValue,
                                           requestedPositions_nm);
    if (canUseCache && cacheMatches(openmmStep,
                                    usesPbc,
                                    boxA_nm,
                                    boxB_nm,
                                    boxC_nm,
                                    lambdaValue,
                                    tabfAlphaValue,
                                    requestedPositions_nm)) {
        biasEnergyKcalPerMol_ = colvars_units::kj_per_mol_to_kcal_per_mol(cachedBiasEnergyKJPerMol_);
        auto& appliedForces = *(modify_atom_applied_forces());
        appliedForces = cachedAppliedForces_;
        return cachedBiasEnergyKJPerMol_;
    }

    int errorCode = updateInput(context, positions);
    errorCode |= cvmodule->calc();
    // OpenMM can re-evaluate forces within a single integrator step.
    // Publish lambda back to the Context only on the first evaluation of that step.
    if (owner_.getAlchemyEnabled() && lastLambdaWritebackStep_ != openmmStep) {
        errorCode |= send_alch_lambda();
        if (errorCode == COLVARS_OK && cvmodule->get_error() == COLVARS_OK) {
            lastLambdaWritebackStep_ = openmmStep;
        }
    }
#ifdef OPENMM_COLVARS_HAS_TABF_ALPHA
    if (owner_.getTabfEnabled() && lastTabfAlphaWritebackStep_ != openmmStep) {
        errorCode |= send_tabf_alpha();
        if (errorCode == COLVARS_OK && cvmodule->get_error() == COLVARS_OK) {
            lastTabfAlphaWritebackStep_ = openmmStep;
        }
    }
#endif
    if (errorCode != COLVARS_OK || cvmodule->get_error() != COLVARS_OK) {
        throw OpenMMException(format_colvars_error(get_error_msgs()));
    }
    const double energyKJPerMol = biasEnergyKJPerMol();
    if (canUseCache) {
        double cacheLambdaValue = lambdaValue;
        if (owner_.getAlchemyEnabled()) {
            cacheLambdaValue = contextParameter(context, owner_.getLambdaParameterName(), lambdaValue);
        }
        double cacheTabfAlphaValue = tabfAlphaValue;
        if (owner_.getTabfEnabled()) {
            cacheTabfAlphaValue = tabfAlpha(context);
        }
        requestedPositions_nm.clear();
        requestedPositions_nm.reserve(atoms_ids.size());
        for (int atomIndex : atoms_ids) {
            if (atomIndex < 0 || atomIndex >= static_cast<int>(positions.size())) {
                invalidateCache();
                return energyKJPerMol;
            }
            requestedPositions_nm.push_back(positions[atomIndex]);
        }
        cachedRequestedPositions_nm_ = std::move(requestedPositions_nm);
        storeCache(openmmStep,
                   usesPbc,
                   boxA_nm,
                   boxB_nm,
                   boxC_nm,
                   cacheLambdaValue,
                   cacheTabfAlphaValue);
    } else {
        invalidateCache();
    }
    return energyKJPerMol;
}

void OpenMMColvarProxy::addForcesToOpenMM(std::vector<Vec3>& forces) const
{
    clearPreviousOpenMMForces(forces);

    const auto& appliedForces = *(get_atom_applied_forces());
    previousForceAtoms_.clear();
    previousForceAtoms_.reserve(atoms_ids.size());
    for (size_t proxyIndex = 0; proxyIndex < atoms_ids.size(); ++proxyIndex) {
        const int openmmAtomIndex = atoms_ids[proxyIndex];
        const cvm::rvector& force_kcal_per_mol_angstrom = appliedForces[proxyIndex];
        constexpr double kForceKcalPerMolAngstromToKJPerMolNm =
            colvars_units::kKJPerKcal * colvars_units::kAngstromsPerNm;
        forces[openmmAtomIndex][0] = force_kcal_per_mol_angstrom.x *
                                     kForceKcalPerMolAngstromToKJPerMolNm;
        forces[openmmAtomIndex][1] = force_kcal_per_mol_angstrom.y *
                                     kForceKcalPerMolAngstromToKJPerMolNm;
        forces[openmmAtomIndex][2] = force_kcal_per_mol_angstrom.z *
                                     kForceKcalPerMolAngstromToKJPerMolNm;
        previousForceAtoms_.push_back(openmmAtomIndex);
    }
}

double OpenMMColvarProxy::biasEnergyKJPerMol() const
{
    return colvars_units::kcal_per_mol_to_kj_per_mol(biasEnergyKcalPerMol_);
}

void OpenMMColvarProxy::updateAtomProperties(int proxyIndex)
{
    const int openmmAtomIndex = atoms_ids[proxyIndex];
    atoms_masses[proxyIndex] = system_.getParticleMass(openmmAtomIndex);
    atoms_charges[proxyIndex] = owner_.getParticleCharge(openmmAtomIndex);
}

void OpenMMColvarProxy::updateAtomIndexMap()
{
    atomIndexById_.clear();
    atomIndexById_.reserve(atoms_ids.size());
    for (size_t index = 0; index < atoms_ids.size(); ++index) {
        atomIndexById_[atoms_ids[index]] = static_cast<int>(index);
    }
}

void OpenMMColvarProxy::updatePbc(ContextImpl& context)
{
    if (!shouldUsePbc(context)) {
        if (boundaries_type != boundaries_non_periodic || pbcLatticeInitialized_) {
            boundaries_type = boundaries_non_periodic;
            reset_pbc_lattice();
            pbcLatticeInitialized_ = false;
        }
        return;
    }

    Vec3 a_nm;
    Vec3 b_nm;
    Vec3 c_nm;
    context.getPeriodicBoxVectors(a_nm, b_nm, c_nm);
    const Boundaries_type newBoundariesType = is_orthorhombic(a_nm, b_nm, c_nm) ? boundaries_pbc_ortho
                                                                                : boundaries_pbc_triclinic;
    if (pbcLatticeInitialized_ && boundaries_type == newBoundariesType &&
        !pbcBoxChanged(a_nm, b_nm, c_nm)) {
        return;
    }

    unit_cell_x.set(colvars_units::nm_to_angstrom(a_nm[0]),
                    colvars_units::nm_to_angstrom(a_nm[1]),
                    colvars_units::nm_to_angstrom(a_nm[2]));
    unit_cell_y.set(colvars_units::nm_to_angstrom(b_nm[0]),
                    colvars_units::nm_to_angstrom(b_nm[1]),
                    colvars_units::nm_to_angstrom(b_nm[2]));
    unit_cell_z.set(colvars_units::nm_to_angstrom(c_nm[0]),
                    colvars_units::nm_to_angstrom(c_nm[1]),
                    colvars_units::nm_to_angstrom(c_nm[2]));

    boundaries_type = newBoundariesType;
    colvarproxy_system::update_pbc_lattice();
    lastPbcA_nm_ = a_nm;
    lastPbcB_nm_ = b_nm;
    lastPbcC_nm_ = c_nm;
    pbcLatticeInitialized_ = true;
}

bool OpenMMColvarProxy::shouldUsePbc(ContextImpl&) const
{
    switch (owner_.getPbcMode()) {
        case ColvarsForce::PbcMode::Enabled:
            return true;
        case ColvarsForce::PbcMode::Disabled:
            return false;
        case ColvarsForce::PbcMode::Auto:
            if (!autoPbcResolved_) {
                cachedAutoPbc_ = resolveAutoPbcMode();
                autoPbcResolved_ = true;
            }
            return cachedAutoPbc_;
    }
    return false;
}

bool OpenMMColvarProxy::resolveAutoPbcMode() const
{
    for (int forceIndex = 0; forceIndex < system_.getNumForces(); ++forceIndex) {
        const Force& force = system_.getForce(forceIndex);
        if (&force == &owner_) {
            continue;
        }
        try {
            if (force.usesPeriodicBoundaryConditions()) {
                return true;
            }
        } catch (const OpenMMException&) {
            continue;
        }
    }
    return false;
}

bool OpenMMColvarProxy::pbcBoxChanged(const Vec3& a_nm, const Vec3& b_nm, const Vec3& c_nm) const
{
    return !same_vec3(a_nm, lastPbcA_nm_) ||
           !same_vec3(b_nm, lastPbcB_nm_) ||
           !same_vec3(c_nm, lastPbcC_nm_);
}

bool OpenMMColvarProxy::buildCacheKey(ContextImpl& context,
                                      const std::vector<Vec3>& positions,
                                      long long& openmmStep,
                                      bool& usesPbc,
                                      Vec3& boxA_nm,
                                      Vec3& boxB_nm,
                                      Vec3& boxC_nm,
                                      double& lambdaValue,
                                      double& tabfAlphaValue,
                                      std::vector<Vec3>& requestedPositions_nm) const
{
    openmmStep = context.getOwner().getStepCount();
    usesPbc = shouldUsePbc(context);
    if (usesPbc) {
        context.getPeriodicBoxVectors(boxA_nm, boxB_nm, boxC_nm);
    } else {
        boxA_nm = Vec3();
        boxB_nm = Vec3();
        boxC_nm = Vec3();
    }

    lambdaValue = 0.0;
    if (owner_.getAlchemyEnabled()) {
        lambdaValue = contextParameter(context, owner_.getLambdaParameterName(), 0.0);
    }
    tabfAlphaValue = owner_.getTabfEnabled() ? tabfAlpha(context) : 1.0;

    requestedPositions_nm.clear();
    requestedPositions_nm.reserve(atoms_ids.size());
    for (int atomIndex : atoms_ids) {
        if (atomIndex < 0 || atomIndex >= static_cast<int>(positions.size())) {
            return false;
        }
        requestedPositions_nm.push_back(positions[atomIndex]);
    }
    return true;
}

bool OpenMMColvarProxy::cacheMatches(long long openmmStep,
                                     bool usesPbc,
                                     const Vec3& boxA_nm,
                                     const Vec3& boxB_nm,
                                     const Vec3& boxC_nm,
                                     double lambdaValue,
                                     double tabfAlphaValue,
                                     const std::vector<Vec3>& requestedPositions_nm) const
{
    if (!cacheValid_ || cachedUsesPbc_ != usesPbc) {
        return false;
    }
    const bool tabfPostStepReadback =
        owner_.getTabfEnabled() && cachedOpenmmStep_ + 1 == openmmStep;
    if (cachedOpenmmStep_ != openmmStep && !tabfPostStepReadback) {
        return false;
    }
    if (!same_scalar(lambdaValue, cachedLambdaValue_)) {
        return false;
    }
    if (!same_scalar(tabfAlphaValue, cachedTabfAlphaValue_)) {
        return false;
    }
    if (cachedAtomIds_ != atoms_ids) {
        return false;
    }
    if (usesPbc &&
        (!same_vec3(cachedPbcA_nm_, boxA_nm) ||
         !same_vec3(cachedPbcB_nm_, boxB_nm) ||
         !same_vec3(cachedPbcC_nm_, boxC_nm))) {
        return false;
    }
    if (requestedPositions_nm.size() != cachedRequestedPositions_nm_.size()) {
        return false;
    }
    for (size_t index = 0; index < requestedPositions_nm.size(); ++index) {
        if (!same_position(requestedPositions_nm[index], cachedRequestedPositions_nm_[index])) {
            return false;
        }
    }
    return cachedAppliedForces_.size() == atoms_ids.size();
}

void OpenMMColvarProxy::storeCache(long long openmmStep,
                                   bool usesPbc,
                                   const Vec3& boxA_nm,
                                   const Vec3& boxB_nm,
                                   const Vec3& boxC_nm,
                                   double lambdaValue,
                                   double tabfAlphaValue)
{
    cachedOpenmmStep_ = openmmStep;
    cachedUsesPbc_ = usesPbc;
    cachedPbcA_nm_ = boxA_nm;
    cachedPbcB_nm_ = boxB_nm;
    cachedPbcC_nm_ = boxC_nm;
    cachedLambdaValue_ = lambdaValue;
    cachedTabfAlphaValue_ = tabfAlphaValue;
    cachedBiasEnergyKJPerMol_ = biasEnergyKJPerMol();
    cachedAtomIds_ = atoms_ids;
    cachedAppliedForces_ = *(get_atom_applied_forces());
    cacheValid_ = true;
}

void OpenMMColvarProxy::invalidateCache()
{
    cacheValid_ = false;
    cachedAtomIds_.clear();
    cachedRequestedPositions_nm_.clear();
    cachedAppliedForces_.clear();
}

void OpenMMColvarProxy::updateColvarsStep(ContextImpl& context)
{
    const long long openmmStep = context.getOwner().getStepCount();
    if (!stepMappingInitialized_) {
        openmmInitialStep_ = openmmStep;
        colvarsInitialStep_ = !owner_.getRestartInputPrefix().empty() || cvmodule->it_restart != 0 ?
            cvmodule->it_restart :
            static_cast<cvm::step_number>(openmmStep);
        cvmodule->it = colvarsInitialStep_;
        stepMappingInitialized_ = true;
        return;
    }

    const long long openmmStepDelta = openmmStep - openmmInitialStep_;
    cvmodule->it = colvarsInitialStep_ + static_cast<cvm::step_number>(openmmStepDelta);
}

void OpenMMColvarProxy::updatePositions(const std::vector<Vec3>& positions)
{
    auto& colvarsPositions = *(modify_atom_positions());
    constexpr double kNmToAngstrom = colvars_units::kAngstromsPerNm;
    for (size_t proxyIndex = 0; proxyIndex < atoms_ids.size(); ++proxyIndex) {
        const Vec3& position_nm = positions[atoms_ids[proxyIndex]];
        colvarsPositions[proxyIndex] = cvm::rvector(
            position_nm[0] * kNmToAngstrom,
            position_nm[1] * kNmToAngstrom,
            position_nm[2] * kNmToAngstrom);
    }
}

void OpenMMColvarProxy::updateTotalForces(ContextImpl& context)
{
    if (!total_force_requested) {
        return;
    }

    context.getForces(systemForceBuffer_);
    if (systemForceBuffer_.size() != static_cast<size_t>(system_.getNumParticles())) {
        set_total_forces_invalid();
        cvmodule->error("Error: OpenMM did not provide total forces required by ABF.\n",
                        COLVARS_ERROR);
        return;
    }

    auto& totalForces = *(modify_atom_total_forces());
    constexpr double kKJPerMolNmToKcalPerMolAngstrom =
        1.0 / (colvars_units::kKJPerKcal * colvars_units::kAngstromsPerNm);
    for (size_t proxyIndex = 0; proxyIndex < atoms_ids.size(); ++proxyIndex) {
        const Vec3& force_kj_per_mol_nm = systemForceBuffer_[atoms_ids[proxyIndex]];
        totalForces[proxyIndex] = cvm::rvector(
            force_kj_per_mol_nm[0] * kKJPerMolNmToKcalPerMolAngstrom,
            force_kj_per_mol_nm[1] * kKJPerMolNmToKcalPerMolAngstrom,
            force_kj_per_mol_nm[2] * kKJPerMolNmToKcalPerMolAngstrom);
    }
    set_total_forces_valid();
}

bool OpenMMColvarProxy::inferTotalForcesSameStep() const
{
    return !config_has_geometrical_eabf(owner_config_text(owner_));
}

void OpenMMColvarProxy::clearPreviousOpenMMForces(std::vector<Vec3>& forces) const
{
    for (int atomIndex : previousForceAtoms_) {
        forces[atomIndex] = Vec3();
    }
}

void OpenMMColvarProxy::captureLambdaDerivative(ContextImpl& context)
{
    currentLambdaDerivativeValid_ = false;
    currentLambdaDerivativeKJPerMol_ = 0.0;
    if (!owner_.getAlchemyEnabled()) {
        return;
    }

    std::map<std::string, double> derivatives;
    context.getEnergyParameterDerivatives(derivatives);
    const auto derivative = derivatives.find(owner_.getLambdaParameterName());
    if (derivative != derivatives.end()) {
        currentLambdaDerivativeKJPerMol_ = derivative->second;
        currentLambdaDerivativeValid_ = true;
        return;
    }

    const auto& parameters = context.getParameters();
    if (parameters.find(owner_.getLambdaDerivativeParameterName()) != parameters.end()) {
        currentLambdaDerivativeKJPerMol_ =
            contextParameter(context, owner_.getLambdaDerivativeParameterName(), 0.0);
        currentLambdaDerivativeValid_ = true;
    }
}

double OpenMMColvarProxy::contextParameter(ContextImpl& context,
                                           const std::string& name,
                                           double defaultValue) const
{
    try {
        return context.getParameter(name);
    } catch (const OpenMMException&) {
        return defaultValue;
    }
}

void OpenMMColvarProxy::setContextParameter(ContextImpl& context,
                                            const std::string& name,
                                            double value) const
{
    try {
        context.setParameter(name, value);
    } catch (const OpenMMException& exception) {
        throw OpenMMException("Colvars-OpenMM failed to set context parameter \"" + name +
                              "\": " + exception.what());
    }
}

double OpenMMColvarProxy::integratorGlobal(ContextImpl& context,
                                           const std::string& name,
                                           double defaultValue) const
{
    auto* customIntegrator = dynamic_cast<CustomIntegrator*>(&context.getIntegrator());
    if (customIntegrator == nullptr) {
        return defaultValue;
    }
    try {
        return customIntegrator->getGlobalVariableByName(name);
    } catch (const OpenMMException&) {
        return defaultValue;
    }
}

bool OpenMMColvarProxy::setIntegratorGlobal(ContextImpl& context,
                                            const std::string& name,
                                            double value) const
{
    auto* customIntegrator = dynamic_cast<CustomIntegrator*>(&context.getIntegrator());
    if (customIntegrator == nullptr) {
        return false;
    }
    try {
        customIntegrator->setGlobalVariableByName(name, value);
        return true;
    } catch (const OpenMMException&) {
        return false;
    }
}

double OpenMMColvarProxy::tabfAlpha(ContextImpl& context) const
{
    return integratorGlobal(context,
                            owner_.getTabfAlphaParameterName(),
                            contextParameter(context, owner_.getTabfAlphaParameterName(), 1.0));
}
