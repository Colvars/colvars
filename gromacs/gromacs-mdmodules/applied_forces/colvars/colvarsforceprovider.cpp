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
 * Implements the force provider for colvars
 *
 * \author Hubert Santuz <hubert.santuz@gmail.com>
 * \ingroup module_applied_forces
 */

#include "colvarsforceprovider.h"

#include <string>

#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/fileio/checkpoint.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdlib/broadcaststructs.h"
#include "gromacs/mdlib/groupcoord.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"


namespace gmx
{

/********************************************************************
 * ColvarsForceProviderState
 */

const std::string ColvarsForceProviderState::nColvarsAtomsName_ = "nColvarsAtoms";

const std::string ColvarsForceProviderState::xOldWholeName_ = "xOldWhole";

const std::string ColvarsForceProviderState::colvarStateFileName_ = "colvarStateFile";

void ColvarsForceProviderState::writeState(KeyValueTreeObjectBuilder kvtBuilder,
                                           const std::string&        identifier) const
{
    writeKvtCheckpointValue(nColvarsAtoms_, nColvarsAtomsName_, identifier, kvtBuilder);

    // Write colvars atoms coords
    auto DoubleArrayAdder = kvtBuilder.addUniformArray<double>(xOldWholeName_);
    for (int i = 0; i < nColvarsAtoms_; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            DoubleArrayAdder.addValue(static_cast<double>(xOldWhole_[i][j]));
        }
    }


    writeKvtCheckpointValue(colvarStateFile_, colvarStateFileName_, identifier, kvtBuilder);
}

void ColvarsForceProviderState::readState(const KeyValueTreeObject& kvtData, const std::string& identifier)
{

    stateRead_ = true;

    readKvtCheckpointValue(compat::make_not_null(&nColvarsAtoms_), nColvarsAtomsName_, identifier, kvtData);


    // Read colvars atoms coords
    auto kvtDoubleArray = kvtData[xOldWholeName_].asArray().values();


    // Make sure the coordinates saved are consistent with the dimensions
    if (kvtDoubleArray.size() % DIM != 0)
    {
        GMX_THROW(InconsistentInputError(
                "Coordinates saved in the checkpoint file are in the wrong format."));
    }

    snew(xOldWhole_, nColvarsAtoms_);
    for (size_t i = 0; i < kvtDoubleArray.size() / DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            xOldWhole_[i][j] = static_cast<real>(kvtDoubleArray[i * DIM + j].cast<double>());
        }
    }

    readKvtCheckpointValue(
            compat::make_not_null(&colvarStateFile_), colvarStateFileName_, identifier, kvtData);
}


/********************************************************************
 * ColvarsForceProvider
 */

ColvarsForceProvider::ColvarsForceProvider(const std::string&       colvarsConfigString,
                                           LocalAtomSetManager*     localAtomSetManager,
                                           PbcType                  pbcType,
                                           double                   simulationTimeStep,
                                           t_atoms                  atoms,
                                           const t_commrec*         cr,
                                           const MDLogger*          logger,
                                           const std::vector<RVec>& colvarsCoords,
                                           const std::string&       outputPrefix,
                                           const std::map<std::string, std::string>& KVTInputs,
                                           const ColvarsForceProviderState&          state,
                                           real                                      ensTemp) :
    ColvarProxyGromacs(colvarsConfigString, atoms, pbcType, logger, MAIN(cr), KVTInputs, ensTemp),
    stateToCheckpoint_(state)
{


    // Total forces on each atom is not available in GROMACS
    total_force_requested = false;

    // Neighbor Search boolean activated during initialization
    gmx_bNS = true;

    // Get GROMACS timestep (picosecond to femtosecond)
    set_integration_timestep(simulationTimeStep * 1000.0);

    output_prefix_str = outputPrefix;


    if (doParsing_)
    {
        colvars->setup_output();
    }


    // MPI initialisation

    // Initialise attributs for the MPI communication
    if (MAIN(cr))
    {
        // Retrieve the number of colvar atoms
        n_colvars_atoms = atoms_ids.size();
    }

    if (PAR(cr))
    {
        // Let the other nodes know the number of colvar atoms and their ids to construct a gmx::LocalAtomSet
        block_bc(cr->mpi_comm_mygroup, n_colvars_atoms);
        atoms_ids.resize(n_colvars_atoms);
        nblock_bc(cr->mpi_comm_mygroup, n_colvars_atoms, atoms_ids.data());

        // Initialise atoms_new_colvar_forces on non-MAIN nodes
        if (!MAIN(cr))
        {
            atoms_new_colvar_forces.resize(n_colvars_atoms);
        }
    }

    // Cast int into Index of the indices for the localAtomSetManager->add() function
    std::vector<Index> index_atoms(atoms_ids.begin(), atoms_ids.end());
    colvars_atoms = std::make_unique<LocalAtomSet>(localAtomSetManager->add(index_atoms));


    snew(x_colvars_unwrapped, n_colvars_atoms);
    snew(xa_shifts, n_colvars_atoms);
    snew(xa_eshifts, n_colvars_atoms);
    snew(f_colvars, n_colvars_atoms);
    snew(xa_old_whole, n_colvars_atoms);


    // Check state status (did we read a cpt file?)
    if (MAIN(cr))
    {
        if (stateToCheckpoint_.stateRead_)
        {
            if (stateToCheckpoint_.nColvarsAtoms_ != n_colvars_atoms)
            {
                cvm::error(
                        "Number of colvars atoms in the .cpt file differs from the one in .tpr "
                        "file");
            }

            // Copy back the last whole positions from the .cpt file
            for (int i = 0; i < n_colvars_atoms; i++)
            {
                copy_rvec(stateToCheckpoint_.xOldWhole_[i], xa_old_whole[i]);
            }

            // Read input state file
            input_buffer_ = stateToCheckpoint_.colvarStateFile_.c_str();
            colvars->setup_input();
        }
        else
        {
            // Initialize state variables
            stateToCheckpoint_.nColvarsAtoms_ = n_colvars_atoms;
            snew(stateToCheckpoint_.xOldWhole_, n_colvars_atoms);

            // Use input coords for the last whole positions.
            for (int i = 0; i < n_colvars_atoms; i++)
            {
                copy_rvec(colvarsCoords[i], xa_old_whole[i]);
            }
        }
    }


    // // Communicate initial coordinates to all processes
    if (PAR(cr))
    {
        nblock_bc(cr->mpi_comm_mygroup, n_colvars_atoms, xa_old_whole);
    }


    if (MAIN(cr) && cvm::debug())
    {
        cvm::log("atoms_ids = " + cvm::to_str(atoms_ids) + "\n");
        cvm::log("atoms_refcount = " + cvm::to_str(atoms_refcount) + "\n");
        cvm::log("positions = " + cvm::to_str(atoms_positions) + "\n");
        cvm::log("total_forces = " + cvm::to_str(atoms_total_forces) + "\n");
        cvm::log("atoms_new_colvar_forces = " + cvm::to_str(atoms_new_colvar_forces) + "\n");
        cvm::log(cvm::line_marker);
        log("done initializing the colvars proxy object.\n");
    }

    if (MAIN(cr))
    {
        cvm::log(cvm::line_marker);
        cvm::log("End colvars Initialization.\n\n");
    }
}

ColvarsForceProvider::~ColvarsForceProvider()
{
    if (doParsing_)
    {
        post_run();
        sfree(stateToCheckpoint_.xOldWhole_);
    }
    sfree(x_colvars_unwrapped);
    sfree(xa_shifts);
    sfree(xa_eshifts);
    sfree(f_colvars);
    sfree(xa_old_whole);
}

void ColvarsForceProvider::calculateForces(const ForceProviderInput& forceProviderInput,
                                           ForceProviderOutput*      forceProviderOutput)
{

    // Construct t_pbc struct
    set_pbc(&gmx_pbc, pbcType_, forceProviderInput.box_);

    const t_commrec* cr = &(forceProviderInput.cr_);
    // Local atom coords
    const gmx::ArrayRef<const gmx::RVec> x = forceProviderInput.x_;
    // Local atom coords (coerced into into old gmx type)
    const rvec* x_pointer = &(x.data()->as_vec());
    const auto& box       = forceProviderInput.box_;

    colvars->it = forceProviderInput.step_;


    // Eventually there needs to be an interface to update local data upon neighbor search
    // We could check if by chance all atoms are in one node, and skip communication
    communicate_group_positions(cr,
                                x_colvars_unwrapped,
                                xa_shifts,
                                xa_eshifts,
                                gmx_bNS,
                                x_pointer,
                                colvars_atoms->numAtomsGlobal(),
                                colvars_atoms->numAtomsLocal(),
                                colvars_atoms->localIndex().data(),
                                colvars_atoms->collectiveIndex().data(),
                                xa_old_whole,
                                box);


    // Communicate_group_positions takes care of removing shifts (unwrapping)
    // in single node jobs, communicate_group_positions() is efficient and adds no overhead

    if (MAIN(cr))
    {
        // On non-MAIN nodes, jump directly to applying the forces

        // backup applied forces if necessary to calculate total forces (if available in future
        // version of Gromacs) if (total_force_requested)
        //  previous_atoms_new_colvar_forces = atoms_new_colvar_forces;

        // Zero the forces on the atoms, so that they can be accumulated by the colvars.
        for (size_t i = 0; i < atoms_new_colvar_forces.size(); i++)
        {
            atoms_new_colvar_forces[i].x         = atoms_new_colvar_forces[i].y =
                    atoms_new_colvar_forces[i].z = 0.0;
        }

        // Get the atom positions from the Gromacs array.
        for (size_t i = 0; i < atoms_ids.size(); i++)
        {
            atoms_positions[i] = cvm::rvector(
                    x_colvars_unwrapped[i][0], x_colvars_unwrapped[i][1], x_colvars_unwrapped[i][2]);
        }

        // // Get total forces if required (if available in future version of Gromacs)
        // if (total_force_requested && cvm::step_relative() > 0) {
        //   for (size_t i = 0; i < atoms_ids.size(); i++) {
        //     size_t aid = atoms_ids[i];
        //     atoms_total_forces[i] = cvm::rvector(f[aid][0], f[aid][1], f[aid][2]);
        //   }
        // }

        bias_energy = 0.0;
        // Call the collective variable module to fill atoms_new_colvar_forces
        if (colvars->calc() != COLVARS_OK)
        {
            cvm::error("Error calling colvars->calc()\n");
        }

        // Copy the forces to C array for broadcasting
        for (int i = 0; i < n_colvars_atoms; i++)
        {
            f_colvars[i][0] = atoms_new_colvar_forces[i].x;
            f_colvars[i][1] = atoms_new_colvar_forces[i].y;
            f_colvars[i][2] = atoms_new_colvar_forces[i].z;
        }

        forceProviderOutput->enerd_.term[F_COM_PULL] += bias_energy;

        // Copy last whole positions into State struct.
        for (int i = 0; i < n_colvars_atoms; i++)
        {
            copy_rvec(xa_old_whole[i], stateToCheckpoint_.xOldWhole_[i]);
        }
    } // MAIN node


    // Broadcast the forces to all the nodes
    if (PAR(cr))
    {
        nblock_bc(cr->mpi_comm_mygroup, n_colvars_atoms, f_colvars);
    }


    const gmx::ArrayRef<gmx::RVec>& f_out = forceProviderOutput->forceWithVirial_.force_;
    matrix                          local_colvars_virial   = { { 0 } };
    const auto&                     localcolvarsIndex      = colvars_atoms->localIndex();
    const auto&                     collectivecolvarsIndex = colvars_atoms->collectiveIndex();
    // Loop through local atoms to aply the colvars forces
    for (gmx::Index l = 0; l < localcolvarsIndex.ssize(); l++)
    {
        /* Get the right index of the local colvars atoms */
        int i_local = localcolvarsIndex[l];
        /* Index of this local atom in the collective colvars atom arrays */
        int i_colvars = collectivecolvarsIndex[l];
        /* Add */
        rvec_inc(f_out[i_local], f_colvars[i_colvars]);
        add_virial_term(local_colvars_virial, f_colvars[i_colvars], x_colvars_unwrapped[i_colvars]);
    }

    forceProviderOutput->forceWithVirial_.addVirialContribution(local_colvars_virial);

    // Re-set the flag for proper update
    gmx_bNS = false;
}

void ColvarsForceProvider::add_virial_term(matrix vir, const rvec& f, const gmx::RVec& x)
{
    for (int j = 0; j < DIM; j++)
    {
        for (int m = 0; m < DIM; m++)
        {
            vir[j][m] -= 0.5 * f[j] * x[m];
        }
    }
}

// Pass restraint energy value for current timestep to MD engine
void ColvarsForceProvider::add_energy(cvm::real energy)
{
    bias_energy += energy;
}


void ColvarsForceProvider::writeCheckpointData(MDModulesWriteCheckpointData checkpointWriting,
                                               const std::string&           moduleName)
{
    colvars->write_restart_string(stateToCheckpoint_.colvarStateFile_);
    stateToCheckpoint_.writeState(checkpointWriting.builder_, moduleName);
}

void ColvarsForceProvider::processAtomsRedistributedSignal(const MDModulesAtomsRedistributedSignal& /*signal*/)
{
    // So far, just update the Neighbor Search boolean for the communicate_group_positions() in calculateForces()
    gmx_bNS = true;
}


} // namespace gmx
