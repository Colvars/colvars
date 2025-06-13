#include "colvargrid_integrate.h"
#include <iostream>
#include <iomanip>
#include <omp.h>


// Helper function to print vector<int>
std::string vec_to_string(const std::vector<int> &vec) {
    std::ostringstream oss;
    oss << "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        oss << vec[i];
        if (i < vec.size() - 1)
            oss << ", ";
    }
    oss << "]";
    return oss.str();
};

colvargrid_integrate::colvargrid_integrate(std::vector<colvar *> &colvars,
                                           std::shared_ptr<colvar_grid_gradient> gradients)
    : colvar_grid_scalar(colvars, gradients, true),
      b_smoothed(false),
      gradients(gradients) {
    // parent class colvar_grid_scalar is constructed with add_extra_bin option set to true
    // hence PMF grid is wider than gradient grid if non-PBC

    if (nd > 1) {
        cvm::main()->cite_feature("Poisson integration of 2D/3D free energy surfaces");
        divergence.resize(nt);

        // Compute inverse of Laplacian diagonal for Jacobi preconditioning
        // For now all code related to preconditioning is commented out
        // until a method better than Jacobi is implemented
        //     cvm::log("Preparing inverse diagonal for preconditioning...\n");
        //     inv_lap_diag.resize(nt);
        //     std::vector<cvm::real> id(nt), lap_col(nt);
        //     for (size_t i = 0; i < nt; i++) {
        //       if (i % (nt / 100) == 0)
        //         cvm::log(cvm::to_str(i));
        //       id[i] = 1.;
        //       atimes(id, lap_col);
        //       id[i] = 0.;
        //       inv_lap_diag[i] = 1. / lap_col[i];
        //     }
        //     cvm::log("Done.\n");
    }
}


colvargrid_integrate::colvargrid_integrate(std::shared_ptr<colvar_grid_gradient> gradients, bool is_weighted)
    : b_smoothed(false), gradients(gradients) {
    nd = gradients->num_variables();
    nx = gradients->number_of_points_vec();
    widths = gradients->widths;
    periodic = gradients->periodic;
    weighted = is_weighted;
    init_computation_nx_nt();
    // Expand grid by 1 bin in non-periodic dimensions
    for (size_t i = 0; i < nd; i++) {
        if (!periodic[i])
            nx[i]++;
        // Shift the grid by half the bin width (values at edges instead of center of bins)
        lower_boundaries.push_back(gradients->lower_boundaries[i].real_value - 0.5 * widths[i]);
    }
    setup(nx);
    if (nd > 1 || weighted) {
        divergence.resize(computation_nt);
    }
    if (weighted) {
        div_border_supplement.resize(computation_nt);
        prepare_divergence_calculation();
        prepare_laplacian_calculation();
    }
    need_to_extrapolate_weighted_solution = false;
    for (size_t i = 0; i < nd; i++) {
        if (!periodic[i]) need_to_extrapolate_weighted_solution = true;
    }
    if (!weighted || !need_to_extrapolate_weighted_solution) {
        computation_grid = this;
    } else {
        computation_grid->periodic = periodic;
        computation_grid->setup(computation_nx);
        // this -> data.clear();
    }
}


int colvargrid_integrate::integrate(const int itmax, const cvm::real &tol, cvm::real &err,
                                    bool verbose) {
    int iter = 0;

    if (nd == 1 && !weighted) {
        cvm::real sum = 0.0;
        cvm::real corr;
        if (periodic[0]) {
            corr = gradients->average(); // Enforce PBC by subtracting average gradient
        } else {
            corr = 0.0;
        }
        std::vector<int> ix;
        // Iterate over valid indices in gradient grid
        for (ix = new_index(); gradients->index_ok(ix); incr(ix)) {
            set_value(ix, sum);
            cvm::real val = gradients->value_output_smoothed(ix, b_smoothed);
            sum += (val - corr) * widths[0];
        }
        if (index_ok(ix)) {
            // This will happen if non-periodic: then PMF grid has one extra bin wrt gradient grid
            // If not, sum should be zero
            set_value(ix, sum);
        }
    } else if (nd <= 3) {
        if (weighted) {
            // // find a nice first starting point
            // divergence.resize(nt);
            // set_div();
            // size_t actual_computation_nt = computation_nt;
            // computation_nt = nt;
            // nr_linbcg_sym(false, divergence, data, 1e-3, 50, iter, err);
            // for (std::vector<int> ix = new_index(); index_ok(ix); incr(ix)) {
            //   bool to_be_considered = true;
            //   std::vector<int> smaller_grid_index = ix;
            //   for (size_t i = 0; i < nd; i++) {
            //     if (ix[i] == 0 or ix[i] == nx[i] - 1) {
            //       to_be_considered = false;
            //       break;
            //     }
            //     else
            //       smaller_grid_index[i] -= 1;
            //   }
            //   if (to_be_considered) {
            //     computation_grid->data[computation_grid->address(smaller_grid_index)] = data[address(ix)];
            //   }
            // }
            divergence.clear();
            // computation_nt = actual_computation_nt;
            // std::cout<< computation_nt << "\n";
            divergence.resize(computation_nt);
            set_weighted_div();
            // TODO: check why it doesn't work when replacing data by computation_grid->data (this was segmentation error)
            // or data or if i use a temporary vector to store this useless laplacian calculation (probably ill-conditionned).
            std::vector<cvm::real> temp = divergence;
            std::vector<cvm::real> temp2(computation_grid->data);
            laplacian_weighted<true>(divergence, temp2);
            temp.clear();
            temp2.clear();
            for (size_t i = 0; i < computation_nt; i++) {
                divergence[i] += div_border_supplement[i];
            }
            div_border_supplement.clear();
        } else {
            set_div();
        }
        optimize_adam(weighted, divergence, computation_grid->data, tol, itmax, iter, err);
        if (verbose)
            cvm::log("Integrated in " + cvm::to_str(iter) + " steps, error: " + cvm::to_str(err));
        if (weighted && need_to_extrapolate_weighted_solution) {
            std::cout << "extrapolating potential" << std::endl;
            // this->data.resize(nt);
            extrapolate_potential();
            std::cout << "potential extrapolated" << std::endl;
        }

        // // DEBUG ###########################
        // auto backup = data;
        // data = divergence;
        // std::ofstream os("div.dat");
        // write_multicol(os);
        // os.close();
        // data = weights;
        // os.open("weights.dat");
        // write_multicol(os);
        // os.close();
        // data = backup;
        // // DEBUG 2 ###########################
        // // Compute terms of the Laplacian matrix
        // std::vector<cvm::real> lap_mat(nt, 0.);

        // std::vector<size_t> cols = {0, 1, 2, 3, 4, 5, nt - 6, nt - 5, nt - 4, nt - 3, nt - 2, nt - 1};

        // for (size_t i = 0; i < cols.size(); i++) {
        //   this->reset();
        //   data[cols[i]] = 1.;
        //   laplacian_weighted<true>(data, lap_mat);
        //   printf("Col  %3li  | ", cols[i]);
        //   for (size_t j = 0; j < cols.size(); j++) {
        //     printf(" %6.1f", lap_mat[cols[j]]);
        //   }
        //   printf("\n");
        // }
        // DEBUG 2 ###########################


        // DEBUG ###########################
        // auto backup = data;
        // data = divergence;
        // std::ofstream os("div.dat");
        // write_multicol(os);
        // os.close();
        // data = weights;
        // os.open("weights.dat");
        // write_multicol(os);
        // os.close();
        // data = backup;
        // DEBUG 2 ###########################
        // Compute terms of the Laplacian matrix
        // std::vector<cvm::real> lap_mat(nt, 0.);

        // std::vector<size_t> cols = { 0, 1, 2, 3, 4, 5, nt-6, nt-5, nt-4, nt-3, nt-2, nt-1 };

        // for (size_t i = 0; i < cols.size(); i++) {
        //   this->reset();
        //   data[cols[i]] = 1.;
        //   laplacian_weighted<true>(data, lap_mat);
        //   printf("Col  %3li  | ", cols[i]);
        //   for (size_t j = 0; j < cols.size(); j++) {
        //     printf(" %6.1f", lap_mat[cols[j]]);
        //   }
        //   printf("\n");
        // }
        // DEBUG 2 ###########################
    } else {
        cvm::error("Cannot integrate PMF in dimension > 3\n");
    }

    return iter;
}


void colvargrid_integrate::set_div() {
    if (nd == 1)
        return;

    for (std::vector<int> ix = new_index(); index_ok(ix); incr(ix)) {
        update_div_local(ix);
    }
}

void colvargrid_integrate::set_weighted_div() {
    sum_count = 0;
    std::vector<int> max_position;
    std::vector<int> min_position;
    sorted_counts = {};

    for (std::vector<int> ix = gradients->new_index(); gradients->index_ok(ix); gradients->incr(ix)) {
        size_t count = gradients->samples->value(ix);
        if (count > 0) {
            insert_into_sorted_list<size_t>(sorted_counts, count);
        }
    }

    upper_threshold_count = sorted_counts[static_cast<int>(sorted_counts.size() * (1 - lambda_max))];
    lower_threshold_count = sorted_counts[static_cast<int>(sorted_counts.size() * lambda_min)];
    sorted_counts.clear();
    size_t n_points = 0;
    for (std::vector<int> ix = gradients->new_index(); gradients->index_ok(ix); gradients->incr(ix)) {
        size_t count = gradients->samples->value(ix);
        if (count < lower_threshold_count) {
            sum_count += lower_threshold_count;
        } else if (count > upper_threshold_count) {
            sum_count += upper_threshold_count;
        } else {
            sum_count += count;
        }
        n_points++;
    }
    m = static_cast<double>(sum_count) / n_points / 1.5;
    std::cout << "mean: " << m * 1.5 << " lower_threshold: " << lower_threshold_count << " upper threshold:" <<
            upper_threshold_count << std::endl;
    regularized_weights.resize(gradients->nt);
    for (std::vector<int> ix = gradients->new_index(); gradients->index_ok(ix);
         gradients->incr(ix)) {
        regularized_weights[gradients->address(ix)] = get_regularized_weight(ix);
    }

    laplacian_coefficients.resize(computation_grid->data.size() * laplacian_stencil.size());
    // laplacian_matrix_test = std::vector<cvm::real>(computation_nt*computation_nt, 0);
    std::vector<int> neighbor_coordinate(nd);
    std::vector<int> reference_point_coordinates(nd);
    std::vector<cvm::real> averaged_normal_vector(nd);
    cvm::real multiplicity = laplacian_stencil.size();
    // Precalculation of the laplacian coefficient

    for (std::vector<int> ix = computation_grid->new_index(); computation_grid->index_ok(ix);
         computation_grid->incr(ix)) {
        for (size_t i = 0; i < laplacian_stencil.size(); i++) {
            std::vector<int> neighbor_relative_position = laplacian_stencil[i];
            for (size_t dim = 0; dim < nd; dim++) {
                neighbor_coordinate[dim] = ix[dim] + neighbor_relative_position[dim];
            }
            cvm::real coefficient = 0;
            //  --> Calling the function makes the code simpler but takes longer to compute
            // coefficient = calculate_weight_sum(neighbor_coordinate,weight_stencil[i]);
            for (std::vector<int> direction: weight_stencil[i]) {
                std::vector<int> weight_coordinate = neighbor_coordinate;
                // Initialize with stencil_point instead of size
                for (size_t n = 0; n < nd && n < direction.size(); n++) {
                    weight_coordinate[n] += direction[n];
                }
                gradients->wrap_detect_edge(weight_coordinate);
                coefficient += regularized_weights[gradients->address(weight_coordinate)] - m;
            }
            coefficient *= weight_counts[i] / pow(2, (nd - 1) * 2);
            coefficient += neighbor_in_classic_laplacian_stencil[i] * m;
            laplacian_coefficients[computation_grid->address(ix) * multiplicity + i] += coefficient;
        }
    }
    for (std::vector<int> ix = computation_grid->new_index(); computation_grid->index_ok(ix);
         computation_grid->incr(ix)) {
        update_weighted_div_local(ix);
    }
}

void colvargrid_integrate::update_div_neighbors(const std::vector<int> &ix0) {
    std::vector<int> ix(ix0);
    int i, j, k;

    // If not periodic, expanded grid ensures that upper neighbors of ix0 are valid grid points
    if (nd == 1) {
        return;
    } else if (nd == 2) {
        update_div_local(ix);
        ix[0]++;
        wrap(ix);
        update_div_local(ix);
        ix[1]++;
        wrap(ix);
        update_div_local(ix);
        ix[0]--;
        wrap(ix);
        update_div_local(ix);
    } else if (nd == 3) {
        for (i = 0; i < 2; i++) {
            ix[1] = ix0[1];
            for (j = 0; j < 2; j++) {
                ix[2] = ix0[2];
                for (k = 0; k < 2; k++) {
                    wrap(ix);
                    update_div_local(ix);
                    ix[2]++;
                }
                ix[1]++;
            }
            ix[0]++;
        }
    }
}


size_t colvargrid_integrate::get_grad(cvm::real *g, std::vector<int> &ix) {
    size_t i;
    bool edge = gradients->wrap_detect_edge(ix); // Detect edge if non-PBC

    if (edge) {
        for (i = 0; i < nd; i++) {
            g[i] = 0.0;
        }
        return 0;
    }

    gradients->vector_value_smoothed(ix, g, b_smoothed);
    if (gradients->samples)
        return gradients->samples->value(ix);
    else
        return 0;
}


size_t colvargrid_integrate::get_grad(std::vector<cvm::real> &g, std::vector<int> &ix) {
    size_t count = gradients->samples->value(ix);
    gradients->vector_value(ix, g);
    return count;
}


void colvargrid_integrate::update_div_local(const std::vector<int> &ix0) {
    const size_t linear_index = address(ix0);
    int i, j, k;
    std::vector<int> ix = ix0;

    if (nd == 2) {
        // gradients at grid points surrounding the current scalar grid point
        cvm::real g00[2], g01[2], g10[2], g11[2];

        get_grad(g11, ix);
        ix[0] = ix0[0] - 1;
        get_grad(g01, ix);
        ix[1] = ix0[1] - 1;
        get_grad(g00, ix);
        ix[0] = ix0[0];
        get_grad(g10, ix);

        divergence[linear_index] = ((g10[0] - g00[0] + g11[0] - g01[0]) / widths[0]
                                    + (g01[1] - g00[1] + g11[1] - g10[1]) / widths[1]) * 0.5;
    } else if (nd == 3) {
        cvm::real gc[24]; // stores 3d gradients in 8 contiguous bins
        int index = 0;

        ix[0] = ix0[0] - 1;
        for (i = 0; i < 2; i++) {
            ix[1] = ix0[1] - 1;
            for (j = 0; j < 2; j++) {
                ix[2] = ix0[2] - 1;
                for (k = 0; k < 2; k++) {
                    get_grad(gc + index, ix);
                    index += 3;
                    ix[2]++;
                }
                ix[1]++;
            }
            ix[0]++;
        }

        divergence[linear_index] =
        ((gc[3 * 4] - gc[0] + gc[3 * 5] - gc[3 * 1] + gc[3 * 6] - gc[3 * 2] + gc[3 * 7] - gc[3 * 3])
         / widths[0]
         + (gc[3 * 2 + 1] - gc[0 + 1] + gc[3 * 3 + 1] - gc[3 * 1 + 1] + gc[3 * 6 + 1] - gc[3 * 4 + 1] + gc[3 * 7 + 1] -
            gc[3 * 5 + 1])
         / widths[1]
         + (gc[3 * 1 + 2] - gc[0 + 2] + gc[3 * 3 + 2] - gc[3 * 2 + 2] + gc[3 * 5 + 2] - gc[3 * 4 + 2] + gc[3 * 7 + 2] -
            gc[3 * 6 + 2])
         / widths[2]) * 0.25;
    }
}


inline size_t min(size_t a, size_t b) { return a < b ? a : b; }

void colvargrid_integrate::prepare_divergence_calculation() {
    surrounding_points_relative_positions.clear();
    size_t n_combinations = pow(2, nd);
    for (size_t i = 0; i < n_combinations; i++) {
        std::string binary = convert_base_two(i, nd);
        std::vector<int> surrounding_point_relative_position = {};
        for (char move_j: binary) {
            surrounding_point_relative_position.push_back(move_j - '0');
        }
        surrounding_points_relative_positions.push_back(surrounding_point_relative_position);
    }
}

void colvargrid_integrate::update_weighted_div_local(const std::vector<int> &ix0)
/*
Updates the divergence at the point ix0
*/
{
    const size_t linear_index = computation_grid->address(ix0);
    cvm::real div_at_point = 0;
    for (std::vector<int> surrounding_point_relative_position:
         surrounding_points_relative_positions) {
        std::vector<int> surrounding_point_coordinates = ix0;
        std::vector<cvm::real> gradient_at_surrounding_point(0, nd);
        for (size_t i = 0; i < nd; i++) {
            surrounding_point_coordinates[i] += surrounding_point_relative_position[i];
        }

        gradients->wrap_detect_edge(surrounding_point_coordinates);
        get_regularized_F(gradient_at_surrounding_point, surrounding_point_coordinates);
        cvm::real weight = regularized_weights[gradients->address(surrounding_point_coordinates)];

        for (size_t i = 0; i < nd; i++) {
            div_at_point +=
                    pow(-1, surrounding_point_relative_position[i] + 1) * gradient_at_surrounding_point[i] * weight /
                    widths[i];
        }
    }
    divergence[linear_index] =
            div_at_point / pow(2, nd - 1);
}

/// Multiplication by sparse matrix representing Laplacian
/// NOTE: Laplacian must be symmetric for solving with CG
void colvargrid_integrate::laplacian(const std::vector<cvm::real> &A, std::vector<cvm::real> &LA) {
    if (nd == 2) {
        // DIMENSION 2

        size_t li, li2;
        int i, j;
        cvm::real fact;
        const cvm::real ffx = 1.0 / (widths[0] * widths[0]);
        const cvm::real ffy = 1.0 / (widths[1] * widths[1]);
        const int h = nx[1];
        const int w = nx[0];
        // offsets for 4 reference points of the Laplacian stencil
        int xm = -h;
        int xp = h;
        int ym = -1;
        int yp = 1;

        // NOTE on performance: this version is slightly sub-optimal because
        // it contains two double loops on the core of the array (for x and y terms)
        // The slightly faster version is in commit 0254cb5a2958cb2e135f268371c4b45fad34866b
        // yet it is much uglier, and probably horrible to extend to dimension 3
        // All terms in the matrix are assigned (=) during the x loops, then updated (+=)
        // with the y (and z) contributions


        // All x components except on x edges
        li = h; // Skip first column

        // Halve the term on y edges (if any) to preserve symmetry of the Laplacian matrix
        // (Long Chen, Finite Difference Methods, UCI, 2017)
        fact = periodic[1] ? 1.0 : 0.5;

        for (i = 1; i < w - 1; i++) {
            // Full range of j, but factor may change on y edges (j == 0 and j == h-1)
            LA[li] = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
            li++;
            for (j = 1; j < h - 1; j++) {
                LA[li] = ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
                li++;
            }
            LA[li] = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
            li++;
        }
        // Edges along x (x components only)
        li = 0L; // Follows left edge
        li2 = h * static_cast<size_t>(w - 1); // Follows right edge
        if (periodic[0]) {
            xm = h * (w - 1);
            xp = h;
            fact = periodic[1] ? 1.0 : 0.5;
            LA[li] = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
            LA[li2] = fact * ffx * (A[li2 - xp] + A[li2 - xm] - 2.0 * A[li2]);
            li++;
            li2++;
            for (j = 1; j < h - 1; j++) {
                LA[li] = ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
                LA[li2] = ffx * (A[li2 - xp] + A[li2 - xm] - 2.0 * A[li2]);
                li++;
                li2++;
            }
            LA[li] = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
            LA[li2] = fact * ffx * (A[li2 - xp] + A[li2 - xm] - 2.0 * A[li2]);
        } else {
            xm = -h;
            xp = h;
            fact = periodic[1] ? 1.0 : 0.5; // Halve in corners in full PBC only
            // lower corner, "j == 0"
            LA[li] = fact * ffx * (A[li + xp] - A[li]);
            LA[li2] = fact * ffx * (A[li2 + xm] - A[li2]);
            li++;
            li2++;
            for (j = 1; j < h - 1; j++) {
                // x gradient (+ y term of laplacian, calculated below)
                LA[li] = ffx * (A[li + xp] - A[li]);
                LA[li2] = ffx * (A[li2 + xm] - A[li2]);
                li++;
                li2++;
            }
            // upper corner, j == h-1
            LA[li] = fact * ffx * (A[li + xp] - A[li]);
            LA[li2] = fact * ffx * (A[li2 + xm] - A[li2]);
        }

        // Now adding all y components
        // All y components except on y edges
        li = 1; // Skip first element (in first row)

        fact = periodic[0] ? 1.0 : 0.5; // for i == 0
        for (i = 0; i < w; i++) {
            // Factor of 1/2 on x edges if non-periodic
            if (i == 1) fact = 1.0;
            if (i == w - 1) fact = periodic[0] ? 1.0 : 0.5;
            for (j = 1; j < h - 1; j++) {
                LA[li] += fact * ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
                li++;
            }
            li += 2; // skip the edges and move to next column
        }
        // Edges along y (y components only)
        li = 0L; // Follows bottom edge
        li2 = h - 1; // Follows top edge
        if (periodic[1]) {
            fact = periodic[0] ? 1.0 : 0.5;
            ym = h - 1;
            yp = 1;
            LA[li] += fact * ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
            LA[li2] += fact * ffy * (A[li2 - yp] + A[li2 - ym] - 2.0 * A[li2]);
            li += h;
            li2 += h;
            for (i = 1; i < w - 1; i++) {
                LA[li] += ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
                LA[li2] += ffy * (A[li2 - yp] + A[li2 - ym] - 2.0 * A[li2]);
                li += h;
                li2 += h;
            }
            LA[li] += fact * ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
            LA[li2] += fact * ffy * (A[li2 - yp] + A[li2 - ym] - 2.0 * A[li2]);
        } else {
            ym = -1;
            yp = 1;
            fact = periodic[0] ? 1.0 : 0.5; // Halve in corners in full PBC only
            // Left corner
            LA[li] += fact * ffy * (A[li + yp] - A[li]);
            LA[li2] += fact * ffy * (A[li2 + ym] - A[li2]);
            li += h;
            li2 += h;
            for (i = 1; i < w - 1; i++) {
                // y gradient (+ x term of laplacian, calculated above)
                LA[li] += ffy * (A[li + yp] - A[li]);
                LA[li2] += ffy * (A[li2 + ym] - A[li2]);
                li += h;
                li2 += h;
            }
            // Right corner
            LA[li] += fact * ffy * (A[li + yp] - A[li]);
            LA[li2] += fact * ffy * (A[li2 + ym] - A[li2]);
        }
    } else if (nd == 3) {
        // DIMENSION 3

        int i, j, k;
        size_t li, li2;
        cvm::real fact = 1.0;
        const cvm::real ffx = 1.0 / (widths[0] * widths[0]);
        const cvm::real ffy = 1.0 / (widths[1] * widths[1]);
        const cvm::real ffz = 1.0 / (widths[2] * widths[2]);
        const int h = nx[2]; // height
        const int d = nx[1]; // depth
        const int w = nx[0]; // width
        // offsets for 6 reference points of the Laplacian stencil
        int xm = -d * h;
        int xp = d * h;
        int ym = -h;
        int yp = h;
        int zm = -1;
        int zp = 1;

        cvm::real factx = periodic[0] ? 1 : 0.5; // factor to be applied on x edges
        cvm::real facty = periodic[1] ? 1 : 0.5; // same for y
        cvm::real factz = periodic[2] ? 1 : 0.5; // same for z
        cvm::real ifactx = 1 / factx;
        cvm::real ifacty = 1 / facty;
        cvm::real ifactz = 1 / factz;

        // All x components except on x edges
        li = d * static_cast<size_t>(h); // Skip left slab
        fact = facty * factz;
        for (i = 1; i < w - 1; i++) {
            for (j = 0; j < d; j++) {
                // full range of y
                if (j == 1) fact *= ifacty;
                if (j == d - 1) fact *= facty;
                LA[li] = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
                li++;
                fact *= ifactz;
                for (k = 1; k < h - 1; k++) {
                    // full range of z
                    LA[li] = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
                    li++;
                }
                fact *= factz;
                LA[li] = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
                li++;
            }
        }
        // Edges along x (x components only)
        li = 0L; // Follows left slab
        li2 = static_cast<size_t>(d) * h * (w - 1); // Follows right slab
        if (periodic[0]) {
            xm = d * h * (w - 1);
            xp = d * h;
            fact = facty * factz;
            for (j = 0; j < d; j++) {
                if (j == 1) fact *= ifacty;
                if (j == d - 1) fact *= facty;
                LA[li] = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
                LA[li2] = fact * ffx * (A[li2 - xp] + A[li2 - xm] - 2.0 * A[li2]);
                li++;
                li2++;
                fact *= ifactz;
                for (k = 1; k < h - 1; k++) {
                    LA[li] = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
                    LA[li2] = fact * ffx * (A[li2 - xp] + A[li2 - xm] - 2.0 * A[li2]);
                    li++;
                    li2++;
                }
                fact *= factz;
                LA[li] = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
                LA[li2] = fact * ffx * (A[li2 - xp] + A[li2 - xm] - 2.0 * A[li2]);
                li++;
                li2++;
            }
        } else {
            xm = -d * h;
            xp = d * h;
            fact = facty * factz;
            for (j = 0; j < d; j++) {
                if (j == 1) fact *= ifacty;
                if (j == d - 1) fact *= facty;
                LA[li] = fact * ffx * (A[li + xp] - A[li]);
                LA[li2] = fact * ffx * (A[li2 + xm] - A[li2]);
                li++;
                li2++;
                fact *= ifactz;
                for (k = 1; k < h - 1; k++) {
                    // x gradient (+ y, z terms of laplacian, calculated below)
                    LA[li] = fact * ffx * (A[li + xp] - A[li]);
                    LA[li2] = fact * ffx * (A[li2 + xm] - A[li2]);
                    li++;
                    li2++;
                }
                fact *= factz;
                LA[li] = fact * ffx * (A[li + xp] - A[li]);
                LA[li2] = fact * ffx * (A[li2 + xm] - A[li2]);
                li++;
                li2++;
            }
        }

        // Now adding all y components
        // All y components except on y edges
        li = h; // Skip first column (in front slab)
        fact = factx * factz;
        for (i = 0; i < w; i++) {
            // full range of x
            if (i == 1) fact *= ifactx;
            if (i == w - 1) fact *= factx;
            for (j = 1; j < d - 1; j++) {
                LA[li] += fact * ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
                li++;
                fact *= ifactz;
                for (k = 1; k < h - 1; k++) {
                    LA[li] += fact * ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
                    li++;
                }
                fact *= factz;
                LA[li] += fact * ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
                li++;
            }
            li += 2 * h; // skip columns in front and back slabs
        }
        // Edges along y (y components only)
        li = 0L; // Follows front slab
        li2 = h * static_cast<size_t>(d - 1); // Follows back slab
        if (periodic[1]) {
            ym = h * (d - 1);
            yp = h;
            fact = factx * factz;
            for (i = 0; i < w; i++) {
                if (i == 1) fact *= ifactx;
                if (i == w - 1) fact *= factx;
                LA[li] += fact * ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
                LA[li2] += fact * ffy * (A[li2 - yp] + A[li2 - ym] - 2.0 * A[li2]);
                li++;
                li2++;
                fact *= ifactz;
                for (k = 1; k < h - 1; k++) {
                    LA[li] += fact * ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
                    LA[li2] += fact * ffy * (A[li2 - yp] + A[li2 - ym] - 2.0 * A[li2]);
                    li++;
                    li2++;
                }
                fact *= factz;
                LA[li] += fact * ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
                LA[li2] += fact * ffy * (A[li2 - yp] + A[li2 - ym] - 2.0 * A[li2]);
                li++;
                li2++;
                li += h * static_cast<size_t>(d - 1);
                li2 += h * static_cast<size_t>(d - 1);
            }
        } else {
            ym = -h;
            yp = h;
            fact = factx * factz;
            for (i = 0; i < w; i++) {
                if (i == 1) fact *= ifactx;
                if (i == w - 1) fact *= factx;
                LA[li] += fact * ffy * (A[li + yp] - A[li]);
                LA[li2] += fact * ffy * (A[li2 + ym] - A[li2]);
                li++;
                li2++;
                fact *= ifactz;
                for (k = 1; k < h - 1; k++) {
                    // y gradient (+ x, z terms of laplacian, calculated above and below)
                    LA[li] += fact * ffy * (A[li + yp] - A[li]);
                    LA[li2] += fact * ffy * (A[li2 + ym] - A[li2]);
                    li++;
                    li2++;
                }
                fact *= factz;
                LA[li] += fact * ffy * (A[li + yp] - A[li]);
                LA[li2] += fact * ffy * (A[li2 + ym] - A[li2]);
                li++;
                li2++;
                li += h * static_cast<size_t>(d - 1);
                li2 += h * static_cast<size_t>(d - 1);
            }
        }

        // Now adding all z components
        // All z components except on z edges
        li = 1; // Skip first element (in bottom slab)
        fact = factx * facty;
        for (i = 0; i < w; i++) {
            // full range of x
            if (i == 1) fact *= ifactx;
            if (i == w - 1) fact *= factx;
            for (k = 1; k < h - 1; k++) {
                LA[li] += fact * ffz * (A[li + zm] + A[li + zp] - 2.0 * A[li]);
                li++;
            }
            fact *= ifacty;
            li += 2; // skip edge slabs
            for (j = 1; j < d - 1; j++) {
                // full range of y
                for (k = 1; k < h - 1; k++) {
                    LA[li] += fact * ffz * (A[li + zm] + A[li + zp] - 2.0 * A[li]);
                    li++;
                }
                li += 2; // skip edge slabs
            }
            fact *= facty;
            for (k = 1; k < h - 1; k++) {
                LA[li] += fact * ffz * (A[li + zm] + A[li + zp] - 2.0 * A[li]);
                li++;
            }
            li += 2; // skip edge slabs
        }
        // Edges along z (z components onlz)
        li = 0; // Follows bottom slab
        li2 = h - 1; // Follows top slab
        if (periodic[2]) {
            zm = h - 1;
            zp = 1;
            fact = factx * facty;
            for (i = 0; i < w; i++) {
                if (i == 1) fact *= ifactx;
                if (i == w - 1) fact *= factx;
                LA[li] += fact * ffz * (A[li + zm] + A[li + zp] - 2.0 * A[li]);
                LA[li2] += fact * ffz * (A[li2 - zp] + A[li2 - zm] - 2.0 * A[li2]);
                li += h;
                li2 += h;
                fact *= ifacty;
                for (j = 1; j < d - 1; j++) {
                    LA[li] += fact * ffz * (A[li + zm] + A[li + zp] - 2.0 * A[li]);
                    LA[li2] += fact * ffz * (A[li2 - zp] + A[li2 - zm] - 2.0 * A[li2]);
                    li += h;
                    li2 += h;
                }
                fact *= facty;
                LA[li] += fact * ffz * (A[li + zm] + A[li + zp] - 2.0 * A[li]);
                LA[li2] += fact * ffz * (A[li2 - zp] + A[li2 - zm] - 2.0 * A[li2]);
                li += h;
                li2 += h;
            }
        } else {
            zm = -1;
            zp = 1;
            fact = factx * facty;
            for (i = 0; i < w; i++) {
                if (i == 1) fact *= ifactx;
                if (i == w - 1) fact *= factx;
                LA[li] += fact * ffz * (A[li + zp] - A[li]);
                LA[li2] += fact * ffz * (A[li2 + zm] - A[li2]);
                li += h;
                li2 += h;
                fact *= ifacty;
                for (j = 1; j < d - 1; j++) {
                    // z gradient (+ x, y terms of laplacian, calculated above)
                    LA[li] += fact * ffz * (A[li + zp] - A[li]);
                    LA[li2] += fact * ffz * (A[li2 + zm] - A[li2]);
                    li += h;
                    li2 += h;
                }
                fact *= facty;
                LA[li] += fact * ffz * (A[li + zp] - A[li]);
                LA[li2] += fact * ffz * (A[li2 + zm] - A[li2]);
                li += h;
                li2 += h;
            }
        }
    }
}


/*
/// Inversion of preconditioner matrix (e.g. diagonal of the Laplacian)
void colvargrid_integrate::asolve(const std::vector<cvm::real> &b, std::vector<cvm::real> &x)
{
  for (size_t i=0; i<int(nt); i++) {
    x[i] = b[i] * inv_lap_diag[i]; // Jacobi preconditioner - little benefit in tests so far
  }
  return;
}*/


/// Multiplication by sparse matrix representing Laplacian
/// NOTE: Laplacian must be symmetric for solving with CG
template<bool initialize_div_supplement>
void colvargrid_integrate::laplacian_weighted(const std::vector<cvm::real> &A, std::vector<cvm::real> &LA) {
    if (initialize_div_supplement) {
        div_border_supplement.resize(divergence.size());
    }
    for (std::vector<int> ix = computation_grid->new_index(); computation_grid->index_ok(ix); computation_grid->
         incr(ix)) {
        LA[computation_grid->address(ix)] = 0;
        if (initialize_div_supplement)
            div_border_supplement[computation_grid->address(ix)] = 0;
         }
    // for testing
    // laplacian_matrix_test = std::vector<cvm::real>(computation_nt*computation_nt, 0);
    if (m_num_threads == 1) {
        for (size_t grid_address = 0; grid_address < computation_grid->nt; grid_address++) {
            std::vector<int> neighbor_coordinate(nd);
            std::vector<int> reference_point_coordinates(nd);
            std::vector<cvm::real> averaged_normal_vector(nd);
            cvm::real multiplicity = laplacian_stencil.size();
            std::vector<int> ix(nd);
            computation_grid->index(grid_address, ix);
            for (size_t i = 0; i < laplacian_stencil.size(); i++) {
                std::vector<int> neighbor_relative_position = laplacian_stencil[i];
                for (size_t dim = 0; dim < nd; dim++) {
                    neighbor_coordinate[dim] = ix[dim] + neighbor_relative_position[dim];
                }
                bool virtual_point = computation_grid->wrap_detect_edge(neighbor_coordinate);
                cvm::real coefficient = laplacian_coefficients[grid_address * multiplicity + i];
                // // Calculation of the laplacian coefficient, when memory is not sufficient to precompute them
                // cvm::real coefficient =0;
                // for (std::vector<int> direction : weight_stencil[i]) {
                //   std::vector<int> weight_coordinate = neighbor_coordinate; // Initialize with stencil_point instead of size
                //   for (size_t n = 0; n < nd && n < direction.size(); n++) {
                //     weight_coordinate[n] += direction[n];
                //   }
                //   gradients->wrap_detect_edge(weight_coordinate);
                //   coefficient += regularized_weights[gradients->address(weight_coordinate)] - m;
                // }
                // coefficient *= weight_counts[i] / pow(2, (nd-1)*2);
                // coefficient+= neighbor_in_classic_laplacian_stencil[i] * m;

                if (!virtual_point) {
                    LA[grid_address] += coefficient * A[computation_grid->address(neighbor_coordinate)];
                    // laplacian_matrix_test[computation_grid->address(ix) * computation_nt + computation_grid->address(neighbor_coordinate)] += coefficient;
                } else {
                    computation_grid->wrap_to_edge(neighbor_coordinate, reference_point_coordinates);
                    LA[grid_address] += coefficient * A[computation_grid->address(reference_point_coordinates)];
                    // laplacian_matrix_test[computation_grid->address(ix) * computation_nt + computation_grid->address(reference_point_coordinates)] += coefficient;

                    if (initialize_div_supplement) {
                        cvm::real div_supplement_term = 0;
                        averaged_normal_vector = compute_averaged_border_normal_gradients(neighbor_coordinate);
                        for (size_t dim = 0; dim < nd; dim++) {
                            div_supplement_term += averaged_normal_vector[dim] * neighbor_relative_position[dim] * widths[
                                dim];
                        }
                        div_border_supplement[grid_address] -= div_supplement_term * coefficient;
                    }
                }
            }
        }
    }
    else{
#ifdef _OPENMP
#pragma omp parallel for // schedule(static, 128)

    for (size_t grid_address = 0; grid_address < computation_grid->nt; grid_address++) {
        std::vector<int> neighbor_coordinate(nd);
        std::vector<int> reference_point_coordinates(nd);
        std::vector<cvm::real> averaged_normal_vector(nd);
        cvm::real multiplicity = laplacian_stencil.size();
        std::vector<int> ix(nd);
        computation_grid->index(grid_address, ix);
        for (size_t i = 0; i < laplacian_stencil.size(); i++) {
            std::vector<int> neighbor_relative_position = laplacian_stencil[i];
            for (size_t dim = 0; dim < nd; dim++) {
                neighbor_coordinate[dim] = ix[dim] + neighbor_relative_position[dim];
            }
            bool virtual_point = computation_grid->wrap_detect_edge(neighbor_coordinate);
            cvm::real coefficient = laplacian_coefficients[grid_address * multiplicity + i];
            // // Calculation of the laplacian coefficient, when memory is not sufficient to precompute them
            // cvm::real coefficient =0;
            // for (std::vector<int> direction : weight_stencil[i]) {
            //   std::vector<int> weight_coordinate = neighbor_coordinate; // Initialize with stencil_point instead of size
            //   for (size_t n = 0; n < nd && n < direction.size(); n++) {
            //     weight_coordinate[n] += direction[n];
            //   }
            //   gradients->wrap_detect_edge(weight_coordinate);
            //   coefficient += regularized_weights[gradients->address(weight_coordinate)] - m;
            // }
            // coefficient *= weight_counts[i] / pow(2, (nd-1)*2);
            // coefficient+= neighbor_in_classic_laplacian_stencil[i] * m;

            if (!virtual_point) {
                LA[grid_address] += coefficient * A[computation_grid->address(neighbor_coordinate)];
                // laplacian_matrix_test[computation_grid->address(ix) * computation_nt + computation_grid->address(neighbor_coordinate)] += coefficient;
            } else {
                computation_grid->wrap_to_edge(neighbor_coordinate, reference_point_coordinates);
                LA[grid_address] += coefficient * A[computation_grid->address(reference_point_coordinates)];
                // laplacian_matrix_test[computation_grid->address(ix) * computation_nt + computation_grid->address(reference_point_coordinates)] += coefficient;

                if (initialize_div_supplement) {
                    cvm::real div_supplement_term = 0;
                    averaged_normal_vector = compute_averaged_border_normal_gradients(neighbor_coordinate);
                    for (size_t dim = 0; dim < nd; dim++) {
                        div_supplement_term += averaged_normal_vector[dim] * neighbor_relative_position[dim] * widths[
                            dim];
                    }
                    div_border_supplement[grid_address] -= div_supplement_term * coefficient;
                }
            }
        }
    }
#else
        cvm::error("multiple threads required in weighted poisson integration, but this binary is not linked "
                   "with a supported threading library.\n");
#endif
}
}

void colvargrid_integrate::prepare_laplacian_calculation() {
    laplacian_stencil.resize(std::pow(3, nd));
    weight_stencil.resize(std::pow(3, nd));
    weight_counts.resize(std::pow(3, nd));
    neighbor_in_classic_laplacian_stencil.resize(std::pow(3, nd));
    for (int i = 0; i < std::pow(3, nd); i++) {
        // for each point in the stencil (for each dimension the relative coordinate can be +1, 0 ,-1)
        std::string base_3 = convert_base_three(i);
        std::vector<int> direction;
        std::vector<std::vector<int> > weights_relative_positions = {
            {}
        }; // relative to the point of the stencil
        double weights_count = 0;
        int dim = 0;
        int number_of_non_zero_coordinates = 0;
        int non_zero_coordinate = -1;
        for (char direction_j: base_3) {
            int displacement_j = direction_j - '0';
            displacement_j -= 1;
            direction.push_back(displacement_j);
            switch (displacement_j) {
                case -1:
                    weights_count += 1.0 / (widths[dim] * widths[dim]);
                    weights_relative_positions =
                            update_weight_relative_positions(weights_relative_positions, std::vector<int>{1});
                    non_zero_coordinate = dim;
                    number_of_non_zero_coordinates++;
                    break;
                case 0:
                    weights_count += -1.0 / (widths[dim] * widths[dim]);
                    weights_relative_positions =
                            update_weight_relative_positions(weights_relative_positions, std::vector<int>{0, 1});
                    break;
                case 1:
                    weights_count += 1.0 / (widths[dim] * widths[dim]);
                    weights_relative_positions =
                            update_weight_relative_positions(weights_relative_positions, std::vector<int>{0});
                    non_zero_coordinate = dim;
                    number_of_non_zero_coordinates++;
                    break;
            }
            dim++;
        }
        // Store computed values in stencil maps
        laplacian_stencil[i] = direction;
        weight_stencil[i] = weights_relative_positions;
        weight_counts[i] = weights_count;

        // Store classic laplacian stencil information
        if (number_of_non_zero_coordinates <= 1) {
            if (non_zero_coordinate != -1)
                neighbor_in_classic_laplacian_stencil[i] =
                        1 / (widths[non_zero_coordinate] * widths[non_zero_coordinate]);
            else {
                float sum = 0;
                for (size_t i = 0; i < nd; i++) {
                    sum -= 2 / (widths[i] * widths[i]);
                }
                neighbor_in_classic_laplacian_stencil[i] = sum;
            }
        } else {
            neighbor_in_classic_laplacian_stencil[i] = 0;
        }
    }
}

void colvargrid_integrate::print_laplacian_preparations() {
    for (int i = 0; i < std::pow(3, nd); i++) {
        std::cout << "Stencil " << i << " is [";
        for (size_t j = 0; j < laplacian_stencil[i].size(); ++j) {
            std::cout << laplacian_stencil[i][j];
            if (j < laplacian_stencil[i].size() - 1)
                std::cout << ", ";
        }
        std::cout << "]" << std::endl;
    }
    std::cout << std::endl;
    std::cout << "weight stencil" << std::endl;
    for (int i = 0; i < std::pow(3, nd); i++) {
        std::cout << "Stencil " << i << " is [";
        for (size_t j = 0; j < weight_stencil[i].size(); ++j) {
            std::cout << vec_to_string(weight_stencil[i][j]);
            if (j < weight_stencil[i].size() - 1)
                std::cout << ", ";
        }
        std::cout << "]" << std::endl;
    }
    std::cout << std::endl;
    std::cout << "weight_counts" << std::endl;
    for (size_t i = 0; i < std::pow(3, nd); i++) {
        std::cout << "Stencil " << i << " is [";
        std::cout << weight_counts[i] << "]" << std::endl;
    }
    std::cout << std::endl;
    std::cout << "neighbor_in_classic_laplacian_stencil" << std::endl;
    for (size_t i = 0; i < std::pow(3, nd); i++) {
        std::cout << "Stencil " << i << " is [";
        std::cout << neighbor_in_classic_laplacian_stencil[i] << "]" << std::endl;
    }
    std::cout << std::endl;
}

std::vector<std::vector<int> > colvargrid_integrate::update_weight_relative_positions(
    std::vector<std::vector<int> > &weights_relative_positions, std::vector<int> direction) {
    std::vector<std::vector<int> > result;

    // For each weight direction and each existing relative position,
    // create a new position by appending the weight direction
    for (int weight_direction: direction) {
        for (const auto &weight_relative_position: weights_relative_positions) {
            // Clone the original position
            std::vector<int> weight_relative_position_clone = weight_relative_position;
            // Append the new direction
            weight_relative_position_clone.push_back(weight_direction);
            // Add to result
            result.push_back(weight_relative_position_clone);
        }
    }
    return result;
}


cvm::real colvargrid_integrate::get_regularized_weight(std::vector<int> &ix) {
    cvm::real regularized_weight;
    size_t count = gradients->samples->value(ix);
    // TODO: put switch here
    if (count < lower_threshold_count) {
        regularized_weight = lower_threshold_count;
    } else if (count > upper_threshold_count) {
        regularized_weight = upper_threshold_count;
    } else {
        regularized_weight = count;
    }
    return regularized_weight;
}

void colvargrid_integrate::get_regularized_F(std::vector<cvm::real> &F, std::vector<int> &ix) {
    // TODO: check if i cannot just use vector_value_smooth_instead/ old get_grad
    F.resize(nd);

    size_t count = get_grad(F, ix);
    float multiplier = 1;
    //TODO: put switch here
    if (count < min_count_F) {
        multiplier = 0;
    } else if (count < max_count_F) {
        multiplier = (count - min_count_F) / (max_count_F - min_count_F);
    }
    for (size_t i = 0; i < nd; i++) {
        F[i] = multiplier * F[i];
    }
}

inline cvm::real colvargrid_integrate::calculate_weight_sum(std::vector<int> stencil_point,
                                                            std::vector<std::vector<int> > directions)
/*
  This function is used to calculate the sum of the weights for a given point of the stencil
  arguments:
    stencil_point: the point of the stencil
    directions: relative positions of the weights
  return:
    the sum of the weights
*/
{
    cvm::real weight_sum = 0;
    for (std::vector<int> direction: directions) {
        std::vector<int> weight_coordinate = stencil_point; // Initialize with stencil_point instead of size
        for (size_t i = 0; i < nd && i < direction.size(); i++) {
            weight_coordinate[i] += direction[i];
        }
        gradients->wrap_detect_edge(weight_coordinate);
        weight_sum += regularized_weights[gradients->address(weight_coordinate)] - m;
    }
    return weight_sum;
}

std::vector<cvm::real> colvargrid_integrate::compute_averaged_border_normal_gradients(
    std::vector<int> virtual_point_coordinates) {
    // bool test = virtual_point_coordinates[0] == 129 && virtual_point_coordinates[1] == 69;
    std::vector<int> reference_point_coordinates(nd, 0); // Initialize with correct size
    gradients->wrap_to_edge(virtual_point_coordinates, reference_point_coordinates);
    std::vector<int> directions_to_average_along;
    std::vector<bool> normal_directions(nd);
    for (size_t i = 0; i < nd; i++) {
        if ((0 <= virtual_point_coordinates[i] && virtual_point_coordinates[i] < computation_nx[i]) || periodic[i]) {
            directions_to_average_along.push_back(i);
            normal_directions[i] = false;
        } else {
            normal_directions[i] = true;
        }
    }
    // Find the positions of the gradients to average
    std::vector<std::vector<int> > gradients_to_average_relative_positions;
    if (directions_to_average_along.size() == 0) {
        std::vector<int> zero_vector(nd, 0);
        gradients_to_average_relative_positions.push_back(zero_vector);
    } else {
        for (int i = 0; i < pow(2, directions_to_average_along.size()); i++) {
            std::vector<int> gradient_to_average_relative_position(nd, 0);
            std::string binary = convert_base_two(i, directions_to_average_along.size());
            for (size_t bit_position = 0; bit_position < directions_to_average_along.size(); bit_position++) {
                gradient_to_average_relative_position[directions_to_average_along[bit_position]] =
                        binary[bit_position] - '0';
            }
            gradients_to_average_relative_positions.push_back(gradient_to_average_relative_position);
        }
    }

    // compute the averaged bordered normal gradient
    std::vector<cvm::real> averaged_bordered_normal_gradient(nd, 0);
    // averaging the gradients
    for (size_t i = 0; i < gradients_to_average_relative_positions.size(); i++) {
        std::vector<int> gradient_position(reference_point_coordinates); // Initialize with reference_point_coordinates
        // if (test){
        //   std::cout<< "gradient_position: " << vec_to_string(gradient_position) << std::endl;
        //   std::cout << vec_to_string(gradients_to_average_relative_positions[i]) << std::endl;
        // }
        for (size_t j = 0; j < nd; j++) {
            gradient_position[j] += gradients_to_average_relative_positions[i][j];
        }
        std::vector<cvm::real> gradient(nd); // Initialize with correct size
        get_regularized_F(gradient, gradient_position);
        for (size_t j = 0; j < nd; j++) {
            averaged_bordered_normal_gradient[j] += gradient[j];
        }
    }
    // only keep the normal directions and average

    for (size_t j = 0; j < nd; j++) {
        if (!normal_directions[j]) {
            averaged_bordered_normal_gradient[j] = 0;
        }
        averaged_bordered_normal_gradient[j] /= gradients_to_average_relative_positions.size();
    }

    return averaged_bordered_normal_gradient;
}

std::string colvargrid_integrate::convert_base_three(int n) {
    std::string result = "";
    // Convert to base 3
    while (n > 0) {
        int remainder = n % 3;
        result.push_back('0' + remainder);
        n /= 3;
    }

    // Handle the case where n is 0
    if (result.empty()) {
        result = "0";
    }

    // Reverse the string (since we built it from right to left)
    reverse(result.begin(), result.end());

    // Pad with leading zeros if necessary
    while (result.size() < nd) {
        result = "0" + result;
    }

    // Truncate if the result has more digits than requested
    if (result.size() > nd) {
        result = result.substr(result.size() - nd);
    }
    return result;
}

std::string colvargrid_integrate::convert_base_two(int n, size_t length) {
    std::string result = "";

    // Convert to base 2
    while (n > 0) {
        int remainder = n % 2;
        result.push_back('0' + remainder);
        n /= 2;
    }

    // Handle the case where n is 0
    if (result.empty()) {
        result = "0";
    }

    // Reverse the string (since we built it from right to left)
    reverse(result.begin(), result.end());

    // Pad with leading zeros if necessary
    while (result.size() < length) {
        result = "0" + result;
    }

    // Truncate if the result has more digits than requested
    if (result.size() > length) {
        result = result.substr(result.size() - length);
    }
    return result;
}


void colvargrid_integrate::nr_linbcg_sym(const bool weighted, const std::vector<cvm::real> &b,
                                         std::vector<cvm::real> &x, const cvm::real &tol,
                                         const int itmax, int &iter, cvm::real &err) {
    cvm::real ak, akden, bk, bkden, bknum, bnrm;
    const cvm::real EPS = 1.0e-14;
    int j;
    std::vector<cvm::real> p(computation_nt), r(computation_nt), z(computation_nt);
    typedef void (colvargrid_integrate::*func_pointer)(const std::vector<double> &,
                                                       std::vector<double> &);
    func_pointer atimes =
            weighted ? &colvargrid_integrate::laplacian_weighted<false> : &colvargrid_integrate::laplacian;

    iter = 0;
    (this->*atimes)(x, r);
    for (j = 0; j < int(computation_nt); j++) {
        r[j] = b[j] - r[j];
    }
    bnrm = l2norm(b);
    if (bnrm < EPS) {
        return; // Target is zero, will break relative error calc
    }
    //   asolve(r,z); // precon
    bkden = 1.0;
    while (iter < itmax) {
        ++iter;
        for (bknum = 0.0, j = 0; j < int(computation_nt); j++) {
            bknum += r[j] * r[j]; // precon: z[j]*r[j]
        }
        if (iter == 1) {
            for (j = 0; j < int(computation_nt); j++) {
                p[j] = r[j]; // precon: p[j] = z[j]
            }
        } else {
            bk = bknum / bkden;
            for (j = 0; j < int(computation_nt); j++) {
                p[j] = bk * p[j] + r[j]; // precon:  bk*p[j] + z[j]
            }
        }
        bkden = bknum;
        (this->*atimes)(p, z);
        for (akden = 0.0, j = 0; j < int(computation_nt); j++) {
            akden += z[j] * p[j];
        }
        ak = bknum / akden;
        for (j = 0; j < int(computation_nt); j++) {
            x[j] += ak * p[j];
            r[j] -= ak * z[j];
        }
        //     asolve(r,z);  // precon
        err = l2norm(r) / bnrm;
        // if (cvm::debug())
        std::cout << "iter=" << std::setw(4) << iter + 1 << std::setw(12) << err << std::endl;
        if (err <= tol)
            break;
    }
}

// --- ADAM Optimizer ---
void colvargrid_integrate::optimize_adam(bool weighted, const std::vector<cvm::real> &b, std::vector<cvm::real> &x,
                                         const cvm::real &tol, const int itmax, int &iter, cvm::real &err
) {
    cvm::real bnrm;
    const cvm::real EPS = 1.0e-14;
    int j;
    std::vector<cvm::real> r(computation_nt, 0);
    x = std::vector<cvm::real>(computation_nt, 0);
    typedef void (colvargrid_integrate::*func_pointer)(const std::vector<double> &,
                                                       std::vector<double> &);
    func_pointer atimes =
            weighted ? &colvargrid_integrate::laplacian_weighted<false> : &colvargrid_integrate::laplacian;

    std::vector<cvm::real> m(computation_nt, 0.0); // First moment vector
    std::vector<cvm::real> v(computation_nt, 0.0); // Second moment vector
    std::vector<cvm::real> v_hat_sqrt_eps(computation_nt);
    cvm::real alpha = 2; // Learning rate
    cvm::real beta1 = 0.9;
    cvm::real beta2 = 0.999;
    cvm::real epsilon = 1e-8;
    cvm::real beta1_t = 1.0; // beta1^t
    cvm::real beta2_t = 1.0; // beta2^t
    cvm::real average_err = 0.0;

    r = std::vector<cvm::real>(computation_nt, 0.0);
    (this->*atimes)(x, r);
    for (j = 0; j < int(computation_nt); j++) {
        // std::cout << x[j] << "  " << r[j] << "  " << b[j] << std::endl;
        r[j] = -r[j] + b[j];
    }
    bnrm = l2norm(b);
    if (bnrm < EPS) {
        return; // Target is zero, will break relative error calc
    }
    //   asolve(r,z); // precon
    bnrm = l2norm(b);
    err = l2norm(r) / bnrm;
    std::cout << "Initial Loss: " << std::fixed << std::setprecision(8) << err << std::endl;

    std::vector<cvm::real> grad(computation_nt, 0.0);
    for (int t = 1; t <= itmax; ++t) {
        (this->*atimes)(x, r);
        for (j = 0; j < int(computation_nt); j++) {
            r[j] = b[j] - r[j];
        }
        //don't forget to update to substract


        beta1_t *= beta1;
        beta2_t *= beta2;

        for (j = 0; j < static_cast<int>(computation_nt); j++) {
            // Update biased first moment estimate
            // m_t = beta1 * m_{t-1} + (1 - beta1) * g_t

            m[j] = beta1 * m[j] + (1 - beta1) * r[j];
            const cvm::real grad_sq_j = r[j] * r[j];
            v[j] = beta2 * v[j] + (1 - beta2) * grad_sq_j;
            const cvm::real m_hat_j = m[j] / (1 - beta1_t);
            const cvm::real v_hat_j = v[j] / (1 - beta2_t);
            v_hat_sqrt_eps[j] = std::sqrt(v_hat_j) + epsilon;
            const cvm::real update_term_j = m_hat_j / v_hat_sqrt_eps[j];
            // if (j==1250)
            //   std::cout << x[j] << "  " << r[j] << "  " << update_term_j << std::endl;
            x[j] -= alpha * update_term_j;
        }
        cvm::real current_err = l2norm(r) / bnrm;
        average_err += current_err;
        // Calculate and print loss
        if (t % 10 == 0 || t == 1) {
            std::cout << "Iteration " << t << ", Loss: " << std::fixed << std::setprecision(8) << current_err <<
                    std::endl;
            if (current_err < tol && t > 10) {
                // t > 10 to avoid early stopping
                std::cout << "Converged " << std::endl;
                break;
            }
            if (average_err / 10.0 > err) {
                alpha *= 0.95;
            }
            err = current_err;
            average_err = 0.0;
        }
    }
}

void colvargrid_integrate::extrapolate_potential() {
    for (std::vector<int> ix = new_index(); index_ok(ix);
         incr(ix)) {
        std::vector<int> corresponding_index_in_small_grid = ix;
        for (size_t i = 0; i < nd; i++) {
            if (!periodic[i]) {
                corresponding_index_in_small_grid[i] = ix[i] - 1;
            }
        }
        std::vector<int> reference_index_in_small_grid(nd, 0);
        bool need_to_extrapolate =
                computation_grid->wrap_to_edge(corresponding_index_in_small_grid, reference_index_in_small_grid);
        cvm::real potential_value = computation_grid->data[computation_grid->address(reference_index_in_small_grid)];
        std::vector<cvm::real> relative_position(nd, 0);

        if (need_to_extrapolate) {
            for (size_t i = 0; i < nd; i++) {
                relative_position[i] = corresponding_index_in_small_grid[i] - reference_index_in_small_grid[i];
            }
            std::vector<cvm::real> averaged_normal_vector = compute_averaged_border_normal_gradients(
                corresponding_index_in_small_grid);
            for (size_t i = 0; i < nd; i++) {
                potential_value += averaged_normal_vector[i] * relative_position[i] * widths[i];
            }
        }
        data[address(ix)] = potential_value;
    }
};

template<typename T>
typename std::vector<T>::iterator colvargrid_integrate::insert_into_sorted_list(
    std::vector<T> &sortedList, const T &value) {
    // Find the first position where the element is not less than value
    auto it = std::lower_bound(sortedList.begin(), sortedList.end(), value);

    // Insert the value at the found position and return iterator to inserted element
    return sortedList.insert(it, value);
}

cvm::real colvargrid_integrate::l2norm(const std::vector<cvm::real> &x) {
    size_t i;
    cvm::real sum = 0.0;
    for (i = 0; i < x.size(); i++)
        sum += x[i] * x[i];
    return sqrt(sum);
}

inline void colvargrid_integrate::reverse(std::string::iterator begin, std::string::iterator end) {
    // Ensure end is after begin
    if (begin >= end) return;

    // Move end iterator back by one since it points one past the last element
    --end;

    // Swap elements from both ends moving toward the center until iterators meet
    while (begin < end) {
        auto temp = *begin;
        *begin = *end;
        *end = temp;
        ++begin;
        --end;
    }
}
