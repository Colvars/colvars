/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Declares the Colvars GROMACS proxy class
 *
 * \author Hubert Santuz <hubert.santuz@gmail.com>
 * \ingroup module_applied_forces
 */

#ifndef GMX_APPLIED_FORCES_COLVARPROXYGROMACS_H
#define GMX_APPLIED_FORCES_COLVARPROXYGROMACS_H

#include "external/colvars/colvaratoms.h"
#include "external/colvars/colvarproxy.h"

#include "gromacs/pbcutil/pbc.h"
#include "gromacs/random/tabulatednormaldistribution.h"
#include "gromacs/random/threefry.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/logger.h"


namespace gmx
{


/*! \internal \brief
 * Implements a GROMACS version of colvarproxy.
 * This class hold for the communication between colvars and GROMACS.
 * 2 child class will inherit from this one: one during pre processing (ColvarsPreProcessor)
 * and one during the simulation (ColvarsForceProvider).
 * Most of the work needed for the communication will be implemented in this class.
 */
class ColvarProxyGromacs : public colvarproxy
{

protected:
    //! Atoms topology
    t_atoms gmx_atoms;

    //! Box infos
    PbcType pbcType_;
    t_pbc   gmx_pbc;

    // GROMACS logger instance
    const MDLogger* logger_ = nullptr;

    //! Activate or not the parsing of the Colvars config file
    bool doParsing_;


    // GROMACS random number generation.
    DefaultRandomEngine           rng; // gromacs random number generator
    TabulatedNormalDistribution<> normal_distribution;


public:
    friend class cvm::atom;

    /*! \brief Construct ColvarProxyGromacs from its parameters
     *
     * \param[in] colvarsConfigString Content of the colvars input file.
     * \param[in] atoms Atoms topology
     * \param[in] pbcType Periodic boundary conditions
     * \param[in] logger GROMACS logger instance
     * \param[in] doParsing Wether the input file should be parsed.
     * \param[in] input_strings Input files stored as string in the KVT
     * \param[in] ensTemp the constant ensemble temperature
     */
    ColvarProxyGromacs(const std::string&                        colvarsConfigString,
                       t_atoms                                   atoms,
                       PbcType                                   pbcType,
                       const MDLogger*                           logger,
                       bool                                      doParsing,
                       const std::map<std::string, std::string>& input_strings,
                       real                                      ensTemp);
    ~ColvarProxyGromacs() override;

    /// Update colvars topology of one atom mass and charge from the GROMACS topology
    void update_atom_properties(int index);

    //! From colvarproxy

    // **************** SYSTEM-WIDE PHYSICAL QUANTITIES ****************
    cvm::real rand_gaussian() override;

    // **************** INPUT/OUTPUT ****************
    /// Print a message to the main log
    void log(std::string const& message) override;
    /// Print a message to the main log and let the rest of the program handle the error
    void error(std::string const& message) override;
    /// Request to set the units used internally by Colvars
    int set_unit_system(std::string const& units_in, bool colvars_defined) override;

    /// Initialize colvars atom from GROMACS topology
    int init_atom(int atom_number) override;

    /*! \brief Check if atom belongs to the global index of atoms
     *  \param[in] atom_number Colvars index of the atom to check
     */
    int check_atom_id(int atom_number) override;

    // **************** PERIODIC BOUNDARY CONDITIONS ****************
    cvm::rvector position_distance(cvm::atom_pos const& pos1, cvm::atom_pos const& pos2) const override;
};

} // namespace gmx

#endif // GMX_APPLIED_FORCES_COLVARPROXYGROMACS_H
