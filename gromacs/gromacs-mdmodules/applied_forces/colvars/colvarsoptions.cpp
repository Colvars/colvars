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
 * Implements options for Colvars.
 */
#include "gmxpre.h"

#include "colvarsoptions.h"

#include <fstream>

#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/optionsection.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/strconvert.h"

#include "colvarspreprocessor.h"


namespace gmx
{

namespace
{

/*! \brief Helper to declare mdp transform rules.
 *
 * Enforces uniform mdp options that are always prepended with the correct
 * string for the densityfitting mdp options.
 *
 * \tparam ToType type to be transformed to
 * \tparam TransformWithFunctionType type of transformation function to be used
 *
 * \param[in] rules KVT transformation rules
 * \param[in] transformationFunction the function to transform the flat kvt tree
 * \param[in] optionTag string tag that describes the mdp option, appended to the
 *                      default string for the density guided simulation
 */
template<class ToType, class TransformWithFunctionType>
void colvarsMdpTransformFromString(IKeyValueTreeTransformRules* rules,
                                   TransformWithFunctionType    transformationFunction,
                                   const std::string&           optionTag)
{
    rules->addRule().from<std::string>("/colvars-" + optionTag).to<ToType>("/colvars/" + optionTag).transformWith(transformationFunction);
}

} // namespace


void ColvarsOptions::initMdpTransform(IKeyValueTreeTransformRules* rules)
{
    const auto& stringIdentityTransform = [](std::string s) { return s; };
    colvarsMdpTransformFromString<bool>(rules, &fromStdString<bool>, c_activeTag_);
    colvarsMdpTransformFromString<std::string>(rules, stringIdentityTransform, colvarsFileNameTag_);
}


void ColvarsOptions::buildMdpOutput(KeyValueTreeObjectBuilder* builder) const
{
    builder->addValue<std::string>("comment-colvars-module", "; Colvars bias");
    builder->addValue<bool>("colvars-" + c_activeTag_, active_);
    builder->addValue<std::string>("comment-colvars-" + colvarsFileNameTag_, "; colvars input file");
    builder->addValue<std::string>("colvars-" + colvarsFileNameTag_, colvarsFileName_);
}


void ColvarsOptions::initMdpOptions(IOptionsContainerWithSections* options)
{
    auto section = options->addSection(OptionSection("colvars"));
    section.addOption(BooleanOption(c_activeTag_.c_str()).store(&active_));
    section.addOption(StringOption(colvarsFileNameTag_.c_str()).store(&colvarsFileName_));
}


void ColvarsOptions::writeInternalParametersToKvt(KeyValueTreeObjectBuilder treeBuilder)
{

    // Copy the content of the colvars input file into a string for latter save in KVT
    std::ifstream f(colvarsFileName_);
    colvarsConfigString =
            std::string((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());

    // Write colvars input file as a string
    treeBuilder.addValue<std::string>("colvars-configString", colvarsConfigString);


    ColvarsPreProcessor colvarsPreProcess(
            colvarsConfigString, gmx_atoms, pbc, logger_, box, x, ensembleTemperature);
    //! Vector with colvars atoms coordinates
    colvarsAtomCoords_ = colvarsPreProcess.getColvarsCoords();

    // Save other colvars input files into the KVT
    if (!colvarsPreProcess.inputStreamsToKVT(treeBuilder, "colvars-inputStreams"))
    {
        GMX_THROW(InternalError("Cannot save colvars input files into the tpr"));
    }

    // Write colvars atoms coords
    auto DoubleArrayAdder = treeBuilder.addUniformArray<double>("colvars-startingCoords");
    for (const auto& indexValue : colvarsAtomCoords_)
    {
        for (int j = 0; j < DIM; j++)
        {
            DoubleArrayAdder.addValue(static_cast<double>(indexValue[j]));
        }
    }

    // Write ensemble temperature
    treeBuilder.addValue<real>("colvars-ensTemp", ensembleTemperature);
}


void ColvarsOptions::readInternalParametersFromKvt(const KeyValueTreeObject& tree)
{

    if (!active_)
    {
        return;
    }


    // Retrieve the content of all inputfiles listed in the KVT as "colvars-inputStreams-filename"
    for (const auto& a : tree.properties())
    {
        std::size_t pos = a.key().find("colvars-inputStreams");
        if (pos != std::string::npos)
        {
            std::string filename = a.key().substr(pos + std::string("colvars-inputStreams").size() + 1);

            inputFiles[filename] = tree[a.key()].cast<std::string>();
        }
    }

    if (!tree.keyExists("colvars-configString"))
    {
        GMX_THROW(InconsistentInputError(
                "Cannot find colvars-configString required for colvars simulation."));
    }
    colvarsConfigString = tree["colvars-configString"].cast<std::string>();


    if (!tree.keyExists("colvars-startingCoords"))
    {
        GMX_THROW(InconsistentInputError(
                "Cannot find colvars-startingCoords required for colvars simulation."));
    }

    auto kvtDoubleArray = tree["colvars-startingCoords"].asArray().values();


    // Make sure the coordinates saved are consistent with the dimensions
    if (kvtDoubleArray.size() % DIM != 0)
    {
        GMX_THROW(InconsistentInputError(
                "Coordinates saved in colvars-startingCoords are in the wrong format."));
    }

    for (size_t i = 0; i < kvtDoubleArray.size() / DIM; i++)
    {
        RVec x;
        for (int j = 0; j < DIM; j++)
        {
            x[j] = static_cast<real>(kvtDoubleArray[i * DIM + j].cast<double>());
        }
        colvarsAtomCoords_.push_back(x);
    }

    ensembleTemperature = tree["colvars-ensTemp"].cast<real>();
}

bool ColvarsOptions::isActive() const
{
    return active_;
}

const std::string& ColvarsOptions::colvarsFileName() const
{
    return colvarsFileName_;
}


void ColvarsOptions::getTopology(gmx_mtop_t* mtop)
{
    gmx_atoms = gmx_mtop_global_atoms(*mtop);
}


const std::string& ColvarsOptions::colvarsInputContent() const
{
    return colvarsConfigString;
}

const std::vector<RVec>& ColvarsOptions::colvarsAtomCoords() const
{
    return colvarsAtomCoords_;
}

const std::string& ColvarsOptions::colvarsOutputPrefix() const
{
    return output_prefix_;
}

const std::map<std::string, std::string>& ColvarsOptions::colvarsInputFiles() const
{
    return inputFiles;
}

const real& ColvarsOptions::colvarsEnsTemp() const
{
    return ensembleTemperature;
}

void ColvarsOptions::processCoordinates(const CoordinatesAndBoxPreprocessed& coord)
{

    x   = coord.coordinates_.unpaddedConstArrayRef();
    pbc = coord.pbc_;
    copy_mat(coord.box_, box);
}

void ColvarsOptions::setLogger(const MDLogger& logger)
{
    logger_ = &logger;
}

void ColvarsOptions::processTprFilename(const MdRunInputFilename& tprFilename)
{
    // Exit if colvars module is not active
    if (!active_)
    {
        return;
    }

    // Provided name should not be empty
    GMX_RELEASE_ASSERT(!tprFilename.mdRunFilename_.empty(),
                       "Filename of the *.tpr simulation file is empty");

    if (!output_prefix_.empty())
    {
        return;
    }

    output_prefix_ = stripExtension(std::filesystem::path(tprFilename.mdRunFilename_).filename());
}

void ColvarsOptions::getTemperature(const EnsembleTemperature& temp)
{
    if (temp.constantEnsembleTemperature_)
    {
        ensembleTemperature = temp.constantEnsembleTemperature_.value();
    }
    else
    {
        ensembleTemperature = -1;
    }
}

} // namespace gmx
