#include <numeric>
#include <algorithm>
// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.


#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarparse.h"
#include "colvar.h"
#include "colvarcomp.h"

colvar::gspath::gspath(std::string const &conf): cvc(conf), atoms(nullptr), reference_frames(0), frame_distances(0) {
    function_type = "gspath";
    // Parse selected atoms
    atoms = parse_group(conf, "atoms");
    // Lookup reference column of PDB
    // Copied from the RMSD class
    std::string reference_column;
    double reference_column_value;
    if (get_keyval(conf, "refPositionsCol", reference_column, std::string(""))) {
        bool found = get_keyval(conf, "refPositionsColValue", reference_column_value, 0.0);
        if (found && reference_column_value == 0.0) {
          cvm::error("Error: refPositionsColValue, "
                     "if provided, must be non-zero.\n");
          return;
        }
    }
    // Lookup all reference frames
    bool has_frames = true;
    size_t num_frames = 1;
    while (has_frames) {
        std::string reference_position_file_lookup = "refPositionFile" + std::to_string(num_frames);
        if (key_lookup(conf, reference_position_file_lookup.c_str())) {
            std::string reference_position_filename;
            get_keyval(conf, reference_position_file_lookup.c_str(), reference_position_filename, std::string(""));
            std::vector<cvm::atom_pos> reference_position(atoms->size());
            cvm::load_coords(reference_position_filename.c_str(), &reference_position, atoms, reference_column, reference_column_value);
            reference_frames.push_back(reference_position);
            ++num_frames;
        } else {
            has_frames = false;
        }
    }
    frame_distances.resize(reference_frames.size());
    frame_index.resize(reference_frames.size());
    std::iota(frame_index.begin(), frame_index.end(), 0);
    v1.resize(atoms->size());
    v2.resize(atoms->size());
    v3.resize(atoms->size());
    dfdv1.resize(atoms->size());
    dfdv2.resize(atoms->size());
    sign = 0;
    M = double(reference_frames.size() - 1);
    m = 1.0;
    // Setup alignment to compute RMSD with respect to reference frames
    for (size_t i_frame = 0; i_frame < reference_frames.size(); ++i_frame) {
        cvm::atom_group* tmp_atoms = parse_group(conf, "atoms");
        // Swipe from the rmsd class
        tmp_atoms->b_center = true;
        tmp_atoms->b_rotate = true;
        tmp_atoms->ref_pos = reference_frames[i_frame];
        tmp_atoms->center_ref_pos();
        tmp_atoms->enable(f_ag_fit_gradients);
        tmp_atoms->rot.request_group1_gradients(tmp_atoms->size());
        tmp_atoms->rot.request_group2_gradients(tmp_atoms->size());
        comp_atoms.push_back(tmp_atoms);
    }
    x.type(colvarvalue::type_scalar);
    // Don't use implicit gradient
    // enable(f_cvc_implicit_gradient);
    // Option for debug gradients
    bool debug_gradients = false;
    get_keyval(conf, "debugGradients", debug_gradients, false);
    if (debug_gradients) enable(f_cvc_debug_gradient);
    get_keyval(conf, "useSecondClosestFrame", use_second_closest_frame, true);
    if (use_second_closest_frame == true) {
        cvm::log(std::string("Geometrical path s(σ) will use the second closest frame to compute s_(m-1)\n"));
    } else {
        cvm::log(std::string("Geometrical path s(σ) will use the neighbouring frame to compute s_(m-1)\n"));
    }
    get_keyval(conf, "useThirdClosestFrame", use_third_closest_frame, false);
    if (use_third_closest_frame == true) {
        cvm::log(std::string("Geometrical path s(σ) will use the third closest frame to compute s_(m+1)\n"));
    } else {
        cvm::log(std::string("Geometrical path s(σ) will use the neighbouring frame to compute s_(m+1)\n"));
    }
    // Logging
    cvm::log(std::string("Geometrical pathCV(s) is initialized.\n"));
    cvm::log(std::string("Geometrical pathCV(s) loaded ") + std::to_string(reference_frames.size()) + std::string(" frames.\n"));
}

void colvar::gspath::update_distances() {
    for (size_t i_frame = 0; i_frame < reference_frames.size(); ++i_frame) {
        cvm::real frame_rmsd = 0.0;
//         comp_atoms[i_frame]->calc_apply_roto_translation();
        for (size_t i_atom = 0; i_atom < atoms->size(); ++i_atom) {
            frame_rmsd += ((*(comp_atoms[i_frame]))[i_atom].pos - reference_frames[i_frame][i_atom]).norm2();
        }
        frame_rmsd /= cvm::real(atoms->size());
        frame_rmsd = cvm::sqrt(frame_rmsd);
        frame_distances[i_frame] = frame_rmsd;
    }
}

void colvar::gspath::calc_value() {
    // Update distance from the frames
    update_distances();
    // Find the closest and the second closest frames
    std::sort(frame_index.begin(), frame_index.end(), [this](size_t i1, size_t i2){return frame_distances[i1] < frame_distances[i2];});
    // Mimic the code in http://www.acmm.nl/ensing/software/PathCV.cpp
    // Determine the sign
    sign = static_cast<long>(frame_index[0]) - static_cast<long>(frame_index[1]);
    if (sign > 1) {
        // sigma(z) is on the left side of the closest frame
        sign = 1;
    } else if (sign < -1) {
        // sigma(z) is on the right side of the closest frame
        sign = -1;
    }
    if (std::abs(static_cast<long>(frame_index[0]) - static_cast<long>(frame_index[1])) > 1) {
        // TODO: Should I throw an error here?
        // NOTE: The math derivation of s and z depends on the assumption that 
        //       the second closest frame is the neibouring frame, which is not 
        //       always true.
        cvm::log(std::string("Warning: Geometrical pathCV(s) relies on the assumption that the second closest frame is the neibouring frame\n"));
        cvm::log(std::string("         Please check your configuration or increase restraint on z(σ)\n"));
        for (size_t i_frame = 0; i_frame < frame_index.size(); ++i_frame) {
            std::string frame_info_line = std::string{"Frame index: "} + std::to_string(frame_index[i_frame]) + std::string{" ; optimal RMSD = "} + std::to_string(frame_distances[frame_index[i_frame]]) + std::string{"\n"};
            cvm::log(frame_info_line);
        }
    }
    min_frame_index_1 = frame_index[0];                                                         // s_m
    min_frame_index_2 = use_second_closest_frame ? frame_index[1] : min_frame_index_1 - sign;   // s_(m-1)
    min_frame_index_3 = use_third_closest_frame ? frame_index[2] : min_frame_index_1 + sign;    // s_(m+1)
    cvm::atom_pos reference_cog_1 = std::accumulate(reference_frames[min_frame_index_1].begin(), reference_frames[min_frame_index_1].end(), cvm::atom_pos(0.0, 0.0, 0.0));
    reference_cog_1 /= reference_frames[min_frame_index_1].size();
    std::vector<cvm::atom_pos> tmp_reference_frame_1(reference_frames[min_frame_index_1].size());
    std::transform(reference_frames[min_frame_index_1].begin(), reference_frames[min_frame_index_1].end(), tmp_reference_frame_1.begin(), [reference_cog_1](const cvm::atom_pos& ai){
        return (ai - reference_cog_1);
    });
    cvm::atom_pos reference_cog_2 = std::accumulate(reference_frames[min_frame_index_2].begin(), reference_frames[min_frame_index_2].end(), cvm::atom_pos(0.0, 0.0, 0.0));
    reference_cog_2 /= reference_frames[min_frame_index_2].size();
    std::vector<cvm::atom_pos> tmp_reference_frame_2(reference_frames[min_frame_index_2].size());
    std::transform(reference_frames[min_frame_index_2].begin(), reference_frames[min_frame_index_2].end(), tmp_reference_frame_2.begin(), [reference_cog_2](const cvm::atom_pos& ai){
        return (ai - reference_cog_2);
    });
    for (size_t i_atom = 0; i_atom < atoms->size(); ++i_atom) {
        // v1 = s_m - z
        v1[i_atom] = reference_frames[min_frame_index_1][i_atom] - (*(comp_atoms[min_frame_index_1]))[i_atom].pos;
        // v2 = z - s_(m-1)
        v2[i_atom] = (*(comp_atoms[min_frame_index_2]))[i_atom].pos - reference_frames[min_frame_index_2][i_atom];
    }
    if (min_frame_index_3 < 0 || min_frame_index_3 > M) {
        // Determine the center of geometry
        rot_v3.calc_optimal_rotation(tmp_reference_frame_1, tmp_reference_frame_2);
        for (size_t i_atom = 0; i_atom < atoms->size(); ++i_atom) {
            v3[i_atom] = rot_v3.q.rotate(tmp_reference_frame_1[i_atom]) - tmp_reference_frame_2[i_atom];
        }
    } else {
        cvm::atom_pos reference_cog_3 = std::accumulate(reference_frames[min_frame_index_3].begin(), reference_frames[min_frame_index_3].end(), cvm::atom_pos(0.0, 0.0, 0.0));
        reference_cog_3 /= reference_frames[min_frame_index_3].size();
        std::vector<cvm::atom_pos> tmp_reference_frame_3(reference_frames[min_frame_index_3].size());
        std::transform(reference_frames[min_frame_index_3].begin(), reference_frames[min_frame_index_3].end(), tmp_reference_frame_3.begin(), [reference_cog_3](const cvm::atom_pos& ai){
            return (ai - reference_cog_3);
        });
        rot_v3.calc_optimal_rotation(tmp_reference_frame_1, tmp_reference_frame_3);
        for (size_t i_atom = 0; i_atom < atoms->size(); ++i_atom) {
            // v3 = s_(m+1) - s_m
            v3[i_atom] = tmp_reference_frame_3[i_atom] - rot_v3.q.rotate(tmp_reference_frame_1[i_atom]);
        }
    }
    // Compute v1v3 and the norm2 of v1, v2 and v3
    v1v3 = 0;
    v1_2 = 0;
    v2_2 = 0;
    v3_2 = 0;
    for (size_t i_atom = 0; i_atom < atoms->size(); ++i_atom) {
        v1v3 += v1[i_atom] * v3[i_atom];
        v1_2 += v1[i_atom] * v1[i_atom];
        v2_2 += v2[i_atom] * v2[i_atom];
        v3_2 += v3[i_atom] * v3[i_atom];
    }
    const cvm::real f = (std::sqrt(v1v3 * v1v3 - v3_2 * (v1_2 - v2_2)) - v1v3) / v3_2;
    // Determine M and m
    m = static_cast<double>(frame_index[0]);
    // Finally compute sigma
    x = m/M + static_cast<double>(sign) * ((f - 1) / (2 * M));
}

void colvar::gspath::calc_gradients() {
    // ATOMS GRADIENTS MUST BE EXPLICIT SET HERE!
    //      Implicit gradients doesn't support fit gradients.
    const cvm::real factor1 = 1.0 / (2.0 * v3_2 * std::sqrt(v1v3 * v1v3 - v3_2 * (v1_2 - v2_2)));
    const cvm::real factor2 = 1.0 / v3_2;
    // Compute the derivative of f with vector v1
    std::transform(v1.begin(), v1.end(), v3.begin(), dfdv1.begin(), [factor1, factor2, this](cvm::atom_pos v1atom, cvm::atom_pos v3atom){
        return (factor1 * (2.0 * v1v3 * v3atom - 2.0 * v3_2 * v1atom) - factor2 * v3atom);
    });
    // Compute the derivative of f with vector v2
    std::transform(v2.begin(), v2.end(), dfdv2.begin(), [factor1, this](cvm::atom_pos v2atom){
        return (factor1 * (2.0 * v3_2 * v2atom));
    });
    cvm::rvector tmp_atom_grad_v1, tmp_atom_grad_v2;
    // dS(v1, v2(r), v3) / dr = ∂S/∂v1 * dv1/dr + ∂S/∂v2 * dv2/dr
    // dv1/dr = [fitting matrix 1][-1, ..., -1]
    // dv2/dr = [fitting matrix 2][1, ..., 1]
    // ∂S/∂v1 = ± (∂f/∂v1) / (2M)
    // ∂S/∂v2 = ± (∂f/∂v2) / (2M)
    // dS(v1, v2(r), v3) / dr = -1.0 * ± (∂f/∂v1) / (2M) + ± (∂f/∂v2) / (2M)
    for (size_t i_atom = 0; i_atom < atoms->size(); ++i_atom) {
        tmp_atom_grad_v1[0] = -1.0 * sign * 0.5 * dfdv1[i_atom][0] / M;
        tmp_atom_grad_v1[1] = -1.0 * sign * 0.5 * dfdv1[i_atom][1] / M;
        tmp_atom_grad_v1[2] = -1.0 * sign * 0.5 * dfdv1[i_atom][2] / M;
        tmp_atom_grad_v2[0] = sign * 0.5 * dfdv2[i_atom][0] / M;
        tmp_atom_grad_v2[1] = sign * 0.5 * dfdv2[i_atom][1] / M;
        tmp_atom_grad_v2[2] = sign * 0.5 * dfdv2[i_atom][2] / M;
        (*(comp_atoms[min_frame_index_1]))[i_atom].grad += tmp_atom_grad_v1;
        (*(comp_atoms[min_frame_index_2]))[i_atom].grad += tmp_atom_grad_v2;
    }
}

void colvar::gspath::apply_force(colvarvalue const &force) {
    // The force applied to this CV is scalar type
    cvm::real const &F = force.real_value;
    (*(comp_atoms[min_frame_index_1])).apply_colvar_force(F);
    (*(comp_atoms[min_frame_index_2])).apply_colvar_force(F);
}

colvar::gzpath::gzpath(std::string const &conf): cvc(conf), atoms(nullptr), reference_frames(0), frame_distances(0) {
    function_type = "gzpath";
    // Parse selected atoms
    atoms = parse_group(conf, "atoms");
    // Lookup reference column of PDB
    // Copied from the RMSD class
    std::string reference_column;
    double reference_column_value;
    if (get_keyval(conf, "refPositionsCol", reference_column, std::string(""))) {
        bool found = get_keyval(conf, "refPositionsColValue", reference_column_value, 0.0);
        if (found && reference_column_value == 0.0) {
          cvm::error("Error: refPositionsColValue, "
                     "if provided, must be non-zero.\n");
          return;
        }
    }
    // Lookup all reference frames
    bool has_frames = true;
    size_t num_frames = 1;
    while (has_frames) {
        std::string reference_position_file_lookup = "refPositionFile" + std::to_string(num_frames);
        if (key_lookup(conf, reference_position_file_lookup.c_str())) {
            std::string reference_position_filename;
            get_keyval(conf, reference_position_file_lookup.c_str(), reference_position_filename, std::string(""));
            std::vector<cvm::atom_pos> reference_position(atoms->size());
            cvm::load_coords(reference_position_filename.c_str(), &reference_position, atoms, reference_column, reference_column_value);
            reference_frames.push_back(reference_position);
            ++num_frames;
        } else {
            has_frames = false;
        }
    }
    frame_distances.resize(reference_frames.size());
    frame_index.resize(reference_frames.size());
    std::iota(frame_index.begin(), frame_index.end(), 0);
    v1.resize(atoms->size());
    v2.resize(atoms->size());
    v3.resize(atoms->size());
    v4.resize(atoms->size());
    dfdv1.resize(atoms->size());
    dfdv2.resize(atoms->size());
    dzdv1.resize(atoms->size());
    dzdv2.resize(atoms->size());
    sign = 0;
    M = double(reference_frames.size() - 1);
    m = 1.0;
    // Setup alignment to compute RMSD with respect to reference frames
    for (size_t i_frame = 0; i_frame < reference_frames.size(); ++i_frame) {
        cvm::atom_group* tmp_atoms = parse_group(conf, "atoms");
        // Swipe from the rmsd class
        tmp_atoms->b_center = true;
        tmp_atoms->b_rotate = true;
        tmp_atoms->ref_pos = reference_frames[i_frame];
        tmp_atoms->center_ref_pos();
        tmp_atoms->enable(f_ag_fit_gradients);
        tmp_atoms->rot.request_group1_gradients(tmp_atoms->size());
        tmp_atoms->rot.request_group2_gradients(tmp_atoms->size());
        comp_atoms.push_back(tmp_atoms);
    }
    x.type(colvarvalue::type_scalar);
    // enable(f_cvc_implicit_gradient);
    // Option for debug gradients
    bool debug_gradients = false;
    get_keyval(conf, "debugGradients", debug_gradients, false);
    if (debug_gradients) enable(f_cvc_debug_gradient);
    get_keyval(conf, "useSecondClosestFrame", use_second_closest_frame, true);
    if (use_second_closest_frame == true) {
        cvm::log(std::string("Geometrical path z(σ) will use the second closest frame to compute s_(m-1)\n"));
    } else {
        cvm::log(std::string("Geometrical path z(σ) will use the neighbouring frame to compute s_(m-1)\n"));
    }
    get_keyval(conf, "useThirdClosestFrame", use_third_closest_frame, false);
    if (use_third_closest_frame == true) {
        cvm::log(std::string("Geometrical path z(σ) will use the third closest frame to compute s_(m+1)\n"));
    } else {
        cvm::log(std::string("Geometrical path z(σ) will use the neighbouring frame to compute s_(m+1)\n"));
    }
    // Logging
    cvm::log(std::string("Geometrical pathCV(z) is initialized.\n"));
    cvm::log(std::string("Geometrical pathCV(z) loaded ") + std::to_string(reference_frames.size()) + std::string(" frames.\n"));
}

void colvar::gzpath::update_distances() {
    for (size_t i_frame = 0; i_frame < reference_frames.size(); ++i_frame) {
        cvm::real frame_rmsd = 0.0;
        for (size_t i_atom = 0; i_atom < atoms->size(); ++i_atom) {
            frame_rmsd += ((*(comp_atoms[i_frame]))[i_atom].pos - reference_frames[i_frame][i_atom]).norm2();
        }
        frame_rmsd /= cvm::real(atoms->size());
        frame_rmsd = cvm::sqrt(frame_rmsd);
        frame_distances[i_frame] = frame_rmsd;
    }
}

void colvar::gzpath::calc_value() {
    // Update distance from the frames
    update_distances();
    // Find the closest and the second closest frames
    std::sort(frame_index.begin(), frame_index.end(), [this](size_t i1, size_t i2){return frame_distances[i1] < frame_distances[i2];});
    // Mimic the code in http://www.acmm.nl/ensing/software/PathCV.cpp
    // Determine the sign
    sign = static_cast<long>(frame_index[0]) - static_cast<long>(frame_index[1]);
    if (sign > 1) {
        // sigma(z) is on the left side of the closest frame
        sign = 1;
    } else if (sign < -1) {
        // sigma(z) is on the right side of the closest frame
        sign = -1;
    }
    if (std::abs(static_cast<long>(frame_index[0]) - static_cast<long>(frame_index[1])) > 1) {
        // TODO: Should I throw an error here?
        // NOTE: The math derivation of s and z depends on the assumption that 
        //       the second closest frame is the neibouring frame, which is not 
        //       always true.
    }
    min_frame_index_1 = frame_index[0];                                                         // s_m
    min_frame_index_2 = use_second_closest_frame ? frame_index[1] : min_frame_index_1 - sign;   // s_(m-1)
    min_frame_index_3 = use_third_closest_frame ? frame_index[2] : min_frame_index_1 + sign;    // s_(m+1)
    cvm::atom_pos reference_cog_1 = std::accumulate(reference_frames[min_frame_index_1].begin(), reference_frames[min_frame_index_1].end(), cvm::atom_pos(0.0, 0.0, 0.0));
    reference_cog_1 /= reference_frames[min_frame_index_1].size();
    std::vector<cvm::atom_pos> tmp_reference_frame_1(reference_frames[min_frame_index_1].size());
    std::transform(reference_frames[min_frame_index_1].begin(), reference_frames[min_frame_index_1].end(), tmp_reference_frame_1.begin(), [reference_cog_1](const cvm::atom_pos& ai){
        return (ai - reference_cog_1);
    });
    cvm::atom_pos reference_cog_2 = std::accumulate(reference_frames[min_frame_index_2].begin(), reference_frames[min_frame_index_2].end(), cvm::atom_pos(0.0, 0.0, 0.0));
    reference_cog_2 /= reference_frames[min_frame_index_2].size();
    std::vector<cvm::atom_pos> tmp_reference_frame_2(reference_frames[min_frame_index_2].size());
    std::transform(reference_frames[min_frame_index_2].begin(), reference_frames[min_frame_index_2].end(), tmp_reference_frame_2.begin(), [reference_cog_2](const cvm::atom_pos& ai){
        return (ai - reference_cog_2);
    });
    rot_v4.calc_optimal_rotation(tmp_reference_frame_1, tmp_reference_frame_2);
    for (size_t i_atom = 0; i_atom < atoms->size(); ++i_atom) {
        v1[i_atom] = reference_frames[min_frame_index_1][i_atom] - (*(comp_atoms[min_frame_index_1]))[i_atom].pos;
        v2[i_atom] = (*(comp_atoms[min_frame_index_2]))[i_atom].pos - reference_frames[min_frame_index_2][i_atom];
        // v4 only computes in gzpath
        // v4 = s_m - s_(m-1)
        v4[i_atom] = rot_v4.q.rotate(tmp_reference_frame_1[i_atom]) - tmp_reference_frame_2[i_atom];
    }
    if (min_frame_index_3 < 0 || min_frame_index_3 > M) {
        // Determine the center of geometry
        rot_v3.calc_optimal_rotation(tmp_reference_frame_1, tmp_reference_frame_2);
        for (size_t i_atom = 0; i_atom < atoms->size(); ++i_atom) {
            v3[i_atom] = rot_v3.q.rotate(tmp_reference_frame_1[i_atom]) - tmp_reference_frame_2[i_atom];
        }
    } else {
        cvm::atom_pos reference_cog_3 = std::accumulate(reference_frames[min_frame_index_3].begin(), reference_frames[min_frame_index_3].end(), cvm::atom_pos(0.0, 0.0, 0.0));
        reference_cog_3 /= reference_frames[min_frame_index_3].size();
        std::vector<cvm::atom_pos> tmp_reference_frame_3(reference_frames[min_frame_index_3].size());
        std::transform(reference_frames[min_frame_index_3].begin(), reference_frames[min_frame_index_3].end(), tmp_reference_frame_3.begin(), [reference_cog_3](const cvm::atom_pos& ai){
            return (ai - reference_cog_3);
        });
        rot_v3.calc_optimal_rotation(tmp_reference_frame_1, tmp_reference_frame_3);
        for (size_t i_atom = 0; i_atom < atoms->size(); ++i_atom) {
            // v3 = s_(m+1) - s_m
            v3[i_atom] = tmp_reference_frame_3[i_atom] - rot_v3.q.rotate(tmp_reference_frame_1[i_atom]);
        }
    }
    // Compute v1v3 and the norm2 of v1, v2 and v3
    v1v3 = 0;
    v1_2 = 0;
    v2_2 = 0;
    v3_2 = 0;
    v4_2 = 0;
    v1v4 = 0;
    for (size_t i_atom = 0; i_atom < atoms->size(); ++i_atom) {
        v1v3 += v1[i_atom] * v3[i_atom];
        v1v4 += v1[i_atom] * v4[i_atom];
        v1_2 += v1[i_atom] * v1[i_atom];
        v2_2 += v2[i_atom] * v2[i_atom];
        v3_2 += v3[i_atom] * v3[i_atom];
        v4_2 += v4[i_atom] * v4[i_atom];
    }
    f = (std::sqrt(v1v3 * v1v3 - v3_2 * (v1_2 - v2_2)) - v1v3) / v3_2;
    dx = 0.5 * (f - 1);
    // z = sqrt((-v1 - dx * v4)^2)
    //   = sqrt(v1^2 + 2dx*(v1⋅v4) + dx * dx * v4^2)
    const cvm::real z_2 = v1_2 + 2 * dx * v1v4 + dx * dx * v4_2;
    z = std::sqrt(z_2);
#ifdef DEBUG_COLVARS_GPATH
    if (cvm::step_absolute() % 1000 == 0) {
        std::cout << "|v1| = " << std::sqrt(v1_2) << " ; |v2| = " << std::sqrt(v2_2) << " ; |v3| = " << std::sqrt(v3_2) << " ; |v4| = " << std::sqrt(v4_2) << '\n';
        for (size_t i_frame = 0; i_frame < frame_index.size(); ++i_frame) {
            std::string frame_info_line = std::string{"Frame index: "} + std::to_string(frame_index[i_frame]) + std::string{" ; optimal RMSD = "} + std::to_string(frame_distances[frame_index[i_frame]]) + std::string{"\n"};
            cvm::log(frame_info_line);
        }
        std::cout << "f = " << f << " ; sign = " << sign << " ; z = " << z << '\n';
    }
#endif
    x = z;
}

void colvar::gzpath::calc_gradients() {
    // The derivative in Ensing's PathCV.cpp seems incomplete
    const cvm::real factor1 = 1.0 / (2.0 * v3_2 * std::sqrt(v1v3 * v1v3 - v3_2 * (v1_2 - v2_2)));
    const cvm::real factor2 = 1.0 / v3_2;
    // Compute the derivative of f with vector v1
    std::transform(v1.begin(), v1.end(), v3.begin(), dfdv1.begin(), [factor1, factor2, this](cvm::atom_pos v1atom, cvm::atom_pos v3atom){
        return (factor1 * (2.0 * v1v3 * v3atom - 2.0 * v3_2 * v1atom) - factor2 * v3atom);
    });
    // Compute the derivative of f with vector v2
    std::transform(v2.begin(), v2.end(), dfdv2.begin(), [factor1, this](cvm::atom_pos v2atom){
        return (factor1 * (2.0 * v3_2 * v2atom));
    });
    // dZ(v1(r), v2(r), v3) / dr = ∂Z/∂v1 * dv1/dr + ∂Z/∂v2 * dv2/dr
    // dv1/dr = [fitting matrix 1][-1, ..., -1]
    // dv2/dr = [fitting matrix 2][1, ..., 1]
    // ∂Z/∂v1 = 1/(2*z) * (2v1 + (f-1)v4 + (v1⋅v4)∂f/∂v1 + v4^2 * 1/4 * 2(f-1) * ∂f/∂v1)
    // ∂Z/∂v2 = 1/(2*z) * ((v1⋅v4)∂f/∂v2 + v4^2 * 1/4 * 2(f-1) * ∂f/∂v2)
    cvm::rvector tmp_atom_grad_v1, tmp_atom_grad_v2;
    for (size_t i_atom = 0; i_atom < atoms->size(); ++i_atom) {
        tmp_atom_grad_v1[0] = -1.0 * (1.0 / (2.0 * z)) * (2.0 * v1[i_atom][0] + (f-1) * v4[i_atom][0] + v1v4 * dfdv1[i_atom][0] + v4_2 * 0.25 * 2.0 * (f-1) * dfdv1[i_atom][0]);
        tmp_atom_grad_v1[1] = -1.0 * (1.0 / (2.0 * z)) * (2.0 * v1[i_atom][1] + (f-1) * v4[i_atom][1] + v1v4 * dfdv1[i_atom][1] + v4_2 * 0.25 * 2.0 * (f-1) * dfdv1[i_atom][1]);
        tmp_atom_grad_v1[2] = -1.0 * (1.0 / (2.0 * z)) * (2.0 * v1[i_atom][2] + (f-1) * v4[i_atom][2] + v1v4 * dfdv1[i_atom][2] + v4_2 * 0.25 * 2.0 * (f-1) * dfdv1[i_atom][2]);
        tmp_atom_grad_v2[0] = (1.0 / (2.0 * z)) * (v1v4 * dfdv2[i_atom][0] + v4_2 * 0.25 * 2.0 * (f-1) * dfdv2[i_atom][0]);
        tmp_atom_grad_v2[1] = (1.0 / (2.0 * z)) * (v1v4 * dfdv2[i_atom][1] + v4_2 * 0.25 * 2.0 * (f-1) * dfdv2[i_atom][1]);
        tmp_atom_grad_v2[2] = (1.0 / (2.0 * z)) * (v1v4 * dfdv2[i_atom][2] + v4_2 * 0.25 * 2.0 * (f-1) * dfdv2[i_atom][2]);
        (*(comp_atoms[min_frame_index_1]))[i_atom].grad += tmp_atom_grad_v1;
        (*(comp_atoms[min_frame_index_2]))[i_atom].grad += tmp_atom_grad_v2;
    }
}

void colvar::gzpath::apply_force(colvarvalue const &force) {
    // The force applied to this CV is scalar type
    cvm::real const &F = force.real_value;
    (*(comp_atoms[min_frame_index_1])).apply_colvar_force(F);
    (*(comp_atoms[min_frame_index_2])).apply_colvar_force(F);
}
