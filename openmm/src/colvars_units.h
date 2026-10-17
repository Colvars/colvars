#ifndef OPENMM_COLVARS_UNITS_H
#define OPENMM_COLVARS_UNITS_H

namespace OpenMM {
namespace colvars_units {

constexpr double kAngstromsPerNm = 10.0;
constexpr double kKJPerKcal = 4.184;

inline double nm_to_angstrom(double position_nm)
{
    return position_nm * kAngstromsPerNm;
}

inline double angstrom_to_nm(double position_angstrom)
{
    return position_angstrom / kAngstromsPerNm;
}

inline double kcal_per_mol_to_kj_per_mol(double energy_kcal_per_mol)
{
    return energy_kcal_per_mol * kKJPerKcal;
}

inline double kj_per_mol_to_kcal_per_mol(double energy_kj_per_mol)
{
    return energy_kj_per_mol / kKJPerKcal;
}

inline double kcal_per_mol_angstrom_to_kj_per_mol_nm(double force_kcal_per_mol_angstrom)
{
    return force_kcal_per_mol_angstrom * kKJPerKcal * kAngstromsPerNm;
}

} // namespace colvars_units
} // namespace OpenMM

#endif
