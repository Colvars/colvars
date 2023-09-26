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
 * Implements the Colvars GROMACS proxy class
 *
 * \author Hubert Santuz <hubert.santuz@gmail.com>
 * \ingroup module_applied_forces
 */

#include "colvarproxygromacs.h"

#include <sstream>


namespace gmx
{

ColvarProxyGromacs::ColvarProxyGromacs(const std::string& colvarsConfigString,
                                       t_atoms            atoms,
                                       PbcType            pbcType,
                                       const MDLogger*    logger,
                                       bool               doParsing,
                                       const std::map<std::string, std::string>& input_strings,
                                       real                                      ensTemp) :
    gmx_atoms(atoms), pbcType_(pbcType), logger_(logger), doParsing_(doParsing)
{

    //! From colvarproxy

    // Retrieve masses and charges from input file
    updated_masses_ = updated_charges_ = true;

    // User-scripted forces are not available in GROMACS
    have_scripts = false;

    angstrom_value_ = 0.1;

    // From Gnu units
    // $ units -ts 'k' 'kJ/mol/K/avogadro'
    // 0.0083144621
    boltzmann_ = 0.0083144621;

    // Get the thermostat temperature.
    set_target_temperature(ensTemp);

    // GROMACS random number generation.
    rng.seed(makeRandomSeed());


    // Read configuration file and set up the proxy during Pre processing
    // and during simulation phase but only on the master node.
    if (doParsing)
    {

        // Retrieve input files stored as string in the KVT
        // Add them to the map of colvars input data.
        for (const auto& [input_name, content] : input_strings)
        {
            input_streams_[input_name] = new std::istringstream(content);
        }

        colvars = new colvarmodule(this);
        cvm::log(cvm::line_marker);
        cvm::log("Start colvars Initialization.");


        colvars->cite_feature("GROMACS engine");
        colvars->cite_feature("Colvars-GROMACS interface");

        if (cvm::debug())
        {
            cvm::log("Initializing the colvars proxy object.\n");
        }

        int error_code = colvarproxy::setup();
        error_code |= colvars->read_config_string(colvarsConfigString);
        error_code |= colvars->update_engine_parameters();
        error_code |= colvars->setup_input();

        if (error_code != COLVARS_OK)
        {
            error("Error when initializing Colvars module.");
        }

        // Citation Reporter
        cvm::log(std::string("\n") + colvars->feature_report(0) + std::string("\n"));

        // TODO: Retrieve step
        // if (step != 0) {
        //     cvm::log("Initializing step number to "+cvm::to_str(step)+".\n");
        // }

        // colvars->it = colvars->it_restart = step;
        colvars->set_initial_step(static_cast<cvm::step_number>(0L));
    }
}


cvm::real ColvarProxyGromacs::rand_gaussian()
{
    return normal_distribution(rng);
}

void ColvarProxyGromacs::log(std::string const& message)
{
    if (logger_)
    {
        std::istringstream is(message);
        std::string        line;
        while (std::getline(is, line))
        {
            GMX_LOG(logger_->info).appendText("colvars: " + line + "\n");
        }
    }
}

void ColvarProxyGromacs::error(std::string const& message)
{
    log(message);
    GMX_THROW(InternalError("Error in collective variables module.\n"));
}


int ColvarProxyGromacs::set_unit_system(std::string const& units_in, bool /*colvars_defined*/)
{
    if (units_in != "gromacs")
    {
        cvm::error(
                "Specified unit system \"" + units_in
                + "\" is unsupported in Gromacs. Supported units are \"gromacs\" (nm, kJ/mol).\n");
        return COLVARS_ERROR;
    }
    return COLVARS_OK;
}


// **************** ATOMS ****************

int ColvarProxyGromacs::check_atom_id(int atom_number)
{
    // GROMACS uses zero-based arrays.
    int const aid = (atom_number - 1);

    if (cvm::debug())
    {
        log("Adding atom " + cvm::to_str(atom_number) + " for collective variables calculation.\n");
    }
    if ((aid < 0) || (aid >= gmx_atoms.nr))
    {
        cvm::error("Error: invalid atom number specified, " + cvm::to_str(atom_number) + "\n",
                   COLVARS_INPUT_ERROR);
        return COLVARS_INPUT_ERROR;
    }

    return aid;
}


int ColvarProxyGromacs::init_atom(int atom_number)
{
    // GROMACS uses zero-based arrays.
    int aid = atom_number - 1;

    for (size_t i = 0; i < atoms_ids.size(); i++)
    {
        if (atoms_ids[i] == aid)
        {
            // this atom id was already recorded
            atoms_refcount[i] += 1;
            return i;
        }
    }

    aid = check_atom_id(atom_number);

    if (aid < 0)
    {
        return COLVARS_INPUT_ERROR;
    }

    int const index = add_atom_slot(aid);
    update_atom_properties(index);
    return index;
}

void ColvarProxyGromacs::update_atom_properties(int index)
{

    // update mass
    double const mass = gmx_atoms.atom[atoms_ids[index]].m;
    if (mass <= 0.001)
    {
        this->log("Warning: near-zero mass for atom " + cvm::to_str(atoms_ids[index] + 1)
                  + "; expect unstable dynamics if you apply forces to it.\n");
    }
    atoms_masses[index] = mass;
    // update charge
    atoms_charges[index] = gmx_atoms.atom[atoms_ids[index]].q;
}

ColvarProxyGromacs::~ColvarProxyGromacs()
{
    if (colvars != nullptr)
    {
        delete colvars;
        colvars = nullptr;
    }
}


cvm::rvector ColvarProxyGromacs::position_distance(cvm::atom_pos const& pos1, cvm::atom_pos const& pos2) const
{
    rvec r1, r2, dr;
    r1[0] = pos1.x;
    r1[1] = pos1.y;
    r1[2] = pos1.z;
    r2[0] = pos2.x;
    r2[1] = pos2.y;
    r2[2] = pos2.z;

    pbc_dx(&gmx_pbc, r2, r1, dr);
    return cvm::atom_pos(dr[0], dr[1], dr[2]);
}


} // namespace gmx
