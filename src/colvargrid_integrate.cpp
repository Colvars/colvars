#include "colvargrid_integrate.h"
#include <iostream>
#include <iomanip>


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
        prepare_divergence_necessary_stencils();
        prepare_laplacian_necessary_stencils();
    }
    need_to_extrapolate_solution = false;
    for (size_t i = 0; i < nd; i++) {
        if (!periodic[i])
            need_to_extrapolate_solution = true;
    }
    if (!need_to_extrapolate_solution) {
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
    } else {
        divergence.clear();
        divergence.resize(computation_nt);
        if (weighted || nd > 3) {
            prepare_calculations();
            // extrapolate_data();
            set_weighted_div();
            std::vector<cvm::real> temp = divergence;
            std::vector<cvm::real> temp2(computation_grid->data);
            linewise_laplacian_weighted =
            precompute ? &colvargrid_integrate::linewise_laplacian_weighted_precomputed<true> :
    &colvargrid_integrate::linewise_laplacian_weighted_otf<true>;
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

        cvm::real average = 0;
        for (size_t i = 0; i < computation_nt; i++) {
            average += divergence[i];
            }
        average = average / static_cast<cvm::real>(computation_nt);
        for (size_t i = 0; i < computation_nt; i++) {
            divergence[i] = divergence[i] - average;
        }

        if (weighted || nd > 3) {
            if (nd > 3) {
                std::cout << "WARNING: Integration of potential of dimension higher than 3 requires a lot of memory." << std::endl;
            }
            nr_linbcg_sym(weighted, divergence, computation_grid->data, tol, itmax, iter, err);
        }
        else
            nr_linbcg_sym(weighted, divergence, computation_grid->data, tol, itmax, iter, err);
        if (verbose)
            cvm::log("Integrated in " + cvm::to_str(iter) + " steps, error: " + cvm::to_str(err));
        if (need_to_extrapolate_solution) {
            std::cout << "extrapolating potential" << std::endl;
            // this->data.resize(nt);
            // if (weighted)
            extrapolate_potential();
            // else
                // extrapolate_potential_unweighted();
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
    }
    // else {
    //     cvm::error("Cannot integrate PMF in dimension > 3\n");
    // }

    return iter;
}


void colvargrid_integrate::set_div() {
    if (nd == 1)
        return;
    for (std::vector<int> ix = computation_grid->new_index(); computation_grid -> index_ok(ix); computation_grid -> incr(ix)) {
        update_div_local(ix);
    }
}

void colvargrid_integrate::set_weighted_div() {
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

size_t colvargrid_integrate::get_grad(std::vector<cvm::real> &g, std::vector<int> ix) {
    size_t count = gradients->samples->value(ix);
    gradients->vector_value(ix, g);
    return count;
}


void colvargrid_integrate::update_div_local(const std::vector<int> &ix0) {
    const size_t linear_index = computation_grid -> address(ix0);
    std::vector<int> ix = ix0;

    if (nd == 2) {
        // gradients at grid points surrounding the current scalar grid point
        std::vector<cvm::real> g00(2,0), g01(2,0), g10(2,0), g11(2,0);


        get_grad(g00, ix);

        ix[0] = ix0[0] + 1;
        computation_grid->wrap_detect_edge(ix);
        get_grad(g10, ix);


        ix[0] = ix0[0] + 1;
        ix[1] = ix0[1] + 1;
        computation_grid->wrap_detect_edge(ix);
        get_grad(g11, ix);


        ix[0] = ix0[0];
        ix[1] = ix0[1] + 1;
        computation_grid->wrap_detect_edge(ix);
        get_grad(g01, ix);

        ix = ix0;

        divergence[linear_index] = ((g10[0] - g00[0] + g11[0] - g01[0]) / widths[0]
                                    + (g01[1] - g00[1] + g11[1] - g10[1]) / widths[1]) * 0.5;
        cvm::real supplement = 0;
        cvm::real supplement_coefficient = 1;
        if ((ix0[0] == 0 || ix0[0] == computation_nx[0] - 1) && !periodic[0]) {
            cvm::real coefficient = (ix[0] == 0) ? 1.0 : -1.0;
            supplement_coefficient = 0.5 * supplement_coefficient;
            supplement += coefficient * (g10[0] + g00[0] + g11[0] + g01[0]) * 0.25 * 2 / widths[0];
        } if ((ix0[1] == 0 || ix0[1] == computation_nx[1] - 1) && !periodic[1]) {
            cvm::real coefficient = (ix[1] == 0) ? 1.0 : -1.0;
            supplement_coefficient = 0.5 * supplement_coefficient;
            supplement += coefficient * (g01[1] + g00[1] + g11[1] + g10[1]) * 0.25 * 2 / widths[1];
        }
        divergence[linear_index] += supplement;
        divergence[linear_index] *= supplement_coefficient;

    } else if (nd == 3) {
        // 8 gradient vectors, each of size 3
        std::vector<cvm::real> g000(3, 0.0), g001(3, 0.0), g010(3, 0.0), g011(3, 0.0);
        std::vector<cvm::real> g100(3, 0.0), g101(3, 0.0), g110(3, 0.0), g111(3, 0.0);


        // Fill the 8 corners, shifting with +1
        get_grad(g000, ix);            // (0,0,0)

        ix[2] = ix0[2] + 1;
        computation_grid->wrap_detect_edge(ix);
        get_grad(g001, ix);            // (0,0,1)

        ix[1] = ix0[1] + 1;
        computation_grid->wrap_detect_edge(ix);
        get_grad(g011, ix);            // (0,1,1)

        ix[2] = ix0[2];
        computation_grid->wrap_detect_edge(ix);
        get_grad(g010, ix);            // (0,1,0)

        ix[0] = ix0[0] + 1;
        ix[1] = ix0[1];
        ix[2] = ix0[2];
        computation_grid->wrap_detect_edge(ix);
        get_grad(g100, ix);            // (1,0,0)

        ix[2] = ix0[2] + 1;
        computation_grid->wrap_detect_edge(ix);
        get_grad(g101, ix);            // (1,0,1)

        ix[1] = ix0[1] + 1;
        computation_grid->wrap_detect_edge(ix);
        get_grad(g111, ix);            // (1,1,1)

        ix[2] = ix0[2];
        computation_grid->wrap_detect_edge(ix);
        get_grad(g110, ix);            // (1,1,0)
        divergence[linear_index] =
        ((g100[0] - g000[0] + g101[0] - g001[0] + g110[0] - g010[0] + g111[0] - g011[0]) / widths[0]
         + (g010[1] - g000[1] + g011[1] - g001[1] + g110[1] - g100[1] + g111[1] - g101[1]) / widths[1]
         + (g001[2] - g000[2] + g011[2] - g010[2] + g101[2] - g100[2] + g111[2] - g110[2]) / widths[2]
        ) * 0.25;
        cvm::real supplement = 0;
        cvm::real supplement_coefficient = 1;

        ix = ix0;
        if ((ix0[0] == 0 || ix0[0] == computation_nx[0] - 1) && !periodic[0]) {
            cvm::real coefficient = (ix[0] == 0) ? 1.0 : -1.0;
            supplement_coefficient *= 0.5;
            supplement = coefficient * (g100[0] + g000[0] + g101[0] + g001[0] + g110[0] + g010[0] + g111[0] + g011[0])
                    * 0.125 * 2 * widths[0];
        } if ((ix0[1] == 0 || ix0[1] == computation_nx[1] - 1) && !periodic[1]) {
            cvm::real coefficient = (ix[1] == 0) ? 1.0 : -1.0;
            supplement_coefficient *= 0.5;
            supplement += coefficient * (g010[1] + g000[1] + g011[1] + g001[1] + g110[1] + g100[1] + g111[1] + g101[1])
                    * widths[1] * 2 * 0.125;
        } if ((ix0[2] == 0 || ix0[2] == computation_nx[2] - 1) && !periodic[2]) {
            cvm::real coefficient = (ix[2] == 0) ? 1.0 : -1.0;
            supplement_coefficient *= 0.5;
            supplement += coefficient * (g001[2] - g000[2] + g011[2] - g010[2] + g101[2] - g100[2] + g111[2] - g110[2])
                    * 0.125 * 2 * widths[2];
        }
        divergence[linear_index] += supplement;
        divergence[linear_index] *= supplement_coefficient;
    }
}

inline size_t min(size_t a, size_t b) { return a < b ? a : b; }

void colvargrid_integrate::prepare_divergence_necessary_stencils() {
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
        std::vector<cvm::real> gradient_at_surrounding_point(nd, 0);
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
        const int h = computation_nx[1];
        const int w = computation_nx[0];
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
        const int h = computation_nx[2]; // height
        const int d = computation_nx[1]; // depth
        const int w = computation_nx[0]; // width
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

template<bool initialize_div_supplement>
void colvargrid_integrate::linewise_laplacian_weighted_precomputed(const std::vector<cvm::real> &A, std::vector<cvm::real> &LA, size_t grid_address) {
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
                bool is_ghost_point = computation_grid->wrap_detect_edge(neighbor_coordinate);
                cvm::real coefficient = laplacian_coefficients[grid_address * multiplicity + i];
                if (!is_ghost_point) {
                    LA[grid_address] += coefficient * A[computation_grid->address(neighbor_coordinate)];
                    // laplacian_matrix_test[computation_grid->address(ix) * computation_nt + computation_grid->address(neighbor_coordinate)] += coefficient;
                } else {
                    computation_grid->wrap_to_edge(neighbor_coordinate, reference_point_coordinates);
                    LA[grid_address] += coefficient * A[computation_grid->address(reference_point_coordinates)];
                    // laplacian_matrix_test[computation_grid->address(ix) * computation_nt + computation_grid->address(reference_point_coordinates)] += coefficient;

                    if (initialize_div_supplement) {
                        cvm::real div_supplement_term = 0;
                        averaged_normal_vector = compute_averaged_border_normal_gradient(neighbor_coordinate);
                        for (size_t dim = 0; dim < nd; dim++) {
                            div_supplement_term += averaged_normal_vector[dim] * neighbor_relative_position[dim] *
                                    widths[dim];
                        }
                        div_border_supplement[grid_address] -= div_supplement_term * coefficient;
                    }
                }
            }
}

template<bool initialize_div_supplement>
void colvargrid_integrate::linewise_laplacian_weighted_otf(const std::vector<cvm::real> &A, std::vector<cvm::real> &LA, size_t grid_address) {
    std::vector<int> neighbor_coordinate(nd);
        std::vector<int> reference_point_coordinates(nd);
        std::vector<cvm::real> averaged_normal_vector(nd);
        std::vector<int> ix(nd);
        computation_grid->index(grid_address, ix);
        cvm::real diagonal_coeff= 0;
        for (size_t i = 0; i < laplacian_stencil.size() - 1; i++) {

            std::vector<int> neighbor_relative_position = laplacian_stencil[i];
            for (size_t dim = 0; dim < nd; dim++) {
                neighbor_coordinate[dim] = ix[dim] + neighbor_relative_position[dim];
            }

            cvm::real coefficient = 0;

            for (std::vector<int> direction: weight_stencil[i]) {
                std::vector<int> weight_coordinate = ix;
                // Initialize with stencil_point instead of size
                for (size_t n = 0; n < nd && n < direction.size(); n++) {
                    weight_coordinate[n] += direction[n];
                }
                gradients->wrap_detect_edge(weight_coordinate);
                coefficient += regularized_weights[gradients->address(weight_coordinate)];
            }
            coefficient *= weight_counts[i];
            diagonal_coeff += coefficient;
            bool is_ghost_point = computation_grid->wrap_detect_edge(neighbor_coordinate);
            if (!is_ghost_point) {
                LA[grid_address] += coefficient * A[computation_grid->address(neighbor_coordinate)];
                // laplacian_matrix_test[computation_grid->address(ix) * computation_nt + computation_grid->address(neighbor_coordinate)] += coefficient;
            } else {
                computation_grid->wrap_to_edge(neighbor_coordinate, reference_point_coordinates);
                LA[grid_address] += coefficient * A[computation_grid->address(reference_point_coordinates)];
                // laplacian_matrix_test[computation_grid->address(ix) * computation_nt + computation_grid->address(reference_point_coordinates)] += coefficient;

                if (initialize_div_supplement) {
                    cvm::real div_supplement_term = 0;
                    averaged_normal_vector = compute_averaged_border_normal_gradient(neighbor_coordinate);
                    for (size_t dim = 0; dim < nd; dim++) {
                        div_supplement_term += averaged_normal_vector[dim] * neighbor_relative_position[dim] *
                                widths[dim];
                    }
                    div_border_supplement[grid_address] -= div_supplement_term * coefficient;
                }
            }
        }
        LA[grid_address] += -diagonal_coeff * A[grid_address];
}

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
    // for testing, uncomment lines involving laplacian_matrix_test, also in linewise_laplacian_weighted
    // laplacian_matrix_test = std::vector<cvm::real>(computation_nt*computation_nt, 0);
    if (m_num_threads == 1) {
        for (size_t grid_address = 0; grid_address < computation_grid->nt; grid_address++) {
            (this->*linewise_laplacian_weighted)(A, LA, grid_address);
        }
    } else {
#ifdef _OPENMP
#pragma omp parallel for // schedule(static, 128)

        for (size_t grid_address = 0; grid_address < computation_grid->nt; grid_address++) {
            (this->*linewise_laplacian_weighted)(A, LA, grid_address);
        }
#else
        cvm::error("multiple threads required in weighted poisson integration, but this binary is not linked "
            "with a supported threading library.\n");
#endif
    }
}

void colvargrid_integrate::prepare_laplacian_necessary_stencils() {
    laplacian_stencil.resize(2 * nd + 1);
    weight_stencil.resize(2 * nd);
    weight_counts.resize(2 * nd);

    for (size_t dim = 0; dim < nd; dim++) {
        std::vector<int> relative_position_left(nd, 0);
        relative_position_left[dim] = -1;
        std::vector<int> relative_position_right(nd, 0);
        relative_position_right[dim] = 1;
        laplacian_stencil[2 * dim] = relative_position_left;
        laplacian_stencil[2 * dim + 1] = relative_position_right;
        weight_counts[2 * dim] = 1 / (widths[dim] * widths[dim]) * (1 / std::pow(2, nd - 1));
        weight_counts[2 * dim + 1] = 1 / (widths[dim] * widths[dim]) * (1 / std::pow(2, nd - 1));
        weight_stencil[2 * dim].resize(std::pow(2, nd - 1));
        weight_stencil[2 * dim + 1].resize(std::pow(2, nd - 1));

        for (int weights_to_average_relative_pos = 0; weights_to_average_relative_pos < std::pow(2, nd - 1);
             weights_to_average_relative_pos++) {
            std::string base_2 = convert_base_two(weights_to_average_relative_pos, nd - 1);
            weight_stencil[2 * dim][weights_to_average_relative_pos].resize(nd);
            weight_stencil[2 * dim + 1][weights_to_average_relative_pos].resize(nd);
            int off_set = 0;
            for (size_t k = 0; k < nd; k++) {
                if (k != dim) {
                    weight_stencil[2 * dim][weights_to_average_relative_pos][k] = base_2[k + off_set] - '0';
                    weight_stencil[2 * dim + 1][weights_to_average_relative_pos][k] = base_2[k + off_set] - '0';
                } else {
                    off_set = -1;
                    weight_stencil[2 * dim][weights_to_average_relative_pos][dim] = 0;
                    weight_stencil[2 * dim + 1][weights_to_average_relative_pos][dim] = 1;
                }
            }
        }
    }
    laplacian_stencil[2 * nd] = std::vector<int>(nd, 0);
    // coefficient for the point where we compute the laplacian
    // is the sum of the stencil's point coefficient. we don't count it in the stencil
}

void colvargrid_integrate::print_laplacian_preparations() {
    for (size_t i = 0; i < laplacian_stencil.size(); i++) {
        std::cout << "Laplacian stencil " << i << " is [";
        for (size_t j = 0; j < laplacian_stencil[i].size(); ++j) {
            std::cout << laplacian_stencil[i][j];
            if (j < laplacian_stencil[i].size() - 1)
                std::cout << ", ";
        }
        std::cout << "]" << std::endl;
    }
    std::cout << std::endl;
    std::cout << "weight stencil" << std::endl;
    for (size_t i = 0; i < weight_stencil.size(); i++) {
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
    for (size_t i = 0; i < weight_counts.size(); i++) {
        std::cout << "Stencil " << i << " is [";
        std::cout << weight_counts[i] << "]" << std::endl;
    }
}

void colvargrid_integrate::prepare_calculations() {
    if (weighted) {
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
        std::cout << " lower_threshold: " << lower_threshold_count << " upper threshold:" <<
                upper_threshold_count << std::endl;
        regularized_weights.resize(gradients->nt);
        for (std::vector<int> ix = gradients->new_index(); gradients->index_ok(ix);
             gradients->incr(ix)) {
            regularized_weights[gradients->address(ix)] = get_regularized_weight(ix);
             }
        std::size_t required = (laplacian_stencil.size() + 2) * computation_nt * sizeof(cvm::real);
        double gigabytes = required / (1024.0 * 1024.0 * 1024.0);
        precompute = gigabytes < 2;
        std::string print_precompute = precompute? "precomputed; precomputing will use " : "computed on the fly, as precomputing would use ";
        std::cout << "Laplacian computation will be " << print_precompute << gigabytes << " GB of memory" << std::endl;
        if (precompute) {
            laplacian_coefficients.clear();
            laplacian_coefficients.resize(computation_grid->data.size() * laplacian_stencil.size(),0);
            // laplacian_matrix_test = std::vector<cvm::real>(computation_nt*computation_nt, 0);
            std::vector<int> neighbor_coordinate(nd);
            std::vector<int> reference_point_coordinates(nd);
            std::vector<cvm::real> averaged_normal_vector(nd);
            cvm::real multiplicity = laplacian_stencil.size();
            // Precalculation of the laplacian coefficient
            cvm::real diagonal_coeff;
            for (std::vector<int> ix = computation_grid->new_index(); computation_grid->index_ok(ix);
                 computation_grid->incr(ix)) {
                diagonal_coeff = 0;
                for (size_t i = 0; i < laplacian_stencil.size() - 1; i++) {
                    std::vector<int> neighbor_relative_position = laplacian_stencil[i];
                    for (size_t dim = 0; dim < nd; dim++) {
                        neighbor_coordinate[dim] = ix[dim] + neighbor_relative_position[dim];
                    }
                    cvm::real coefficient = 0;
                    //  --> Calling the function makes the code simpler but takes longer to compute
                    // coefficient = calculate_weight_sum(ix,weight_stencil[i]);
                    for (std::vector<int> direction: weight_stencil[i]) {
                        std::vector<int> weight_coordinate = ix;
                        // Initialize with stencil_point instead of size
                        for (size_t n = 0; n < nd && n < direction.size(); n++) {
                            weight_coordinate[n] += direction[n];
                        }
                        gradients->wrap_detect_edge(weight_coordinate);
                        coefficient += regularized_weights[gradients->address(weight_coordinate)];
                    }
                    coefficient *= weight_counts[i];
                    laplacian_coefficients[computation_grid->address(ix) * multiplicity + i] += coefficient;
                    diagonal_coeff += coefficient;
                }
                laplacian_coefficients[computation_grid->address(ix) * multiplicity + (2 * nd)] += -diagonal_coeff;
            }
        }
    }
}

cvm::real colvargrid_integrate::get_regularized_weight(std::vector<int> &ix) {
    cvm::real regularized_weight;
    size_t count = gradients->samples->value(ix);
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
    // F.resize(nd);

    size_t count = get_grad(F, ix);
    float multiplier = 1;
    if (count < min_count_F) {
        multiplier = 0;
    } else if (count < max_count_F) {
        multiplier = (count - min_count_F) / (max_count_F - min_count_F);
    }
    for (size_t i = 0; i < nd; i++) {
        F[i] = multiplier * F[i];
    }
}

std::vector<cvm::real> colvargrid_integrate::compute_averaged_border_normal_gradient(
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
    linewise_laplacian_weighted =
            precompute ? &colvargrid_integrate::linewise_laplacian_weighted_precomputed<false> :
    &colvargrid_integrate::linewise_laplacian_weighted_otf<false>;
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
    linewise_laplacian_weighted =
            precompute ? &colvargrid_integrate::linewise_laplacian_weighted_precomputed<false> :
    &colvargrid_integrate::linewise_laplacian_weighted_otf<false>;
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
    for (iter = 1; iter <= itmax; ++iter) {
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
        if (iter % 10 == 0 || iter == 1) {
            std::cout << "Iteration " << iter << ", Loss: " << std::fixed << std::setprecision(8) << current_err <<
                    std::endl;
            if (current_err < tol && iter > 10) {
                // t > 10 to avoid early stopping
                std::cout << "Converged " << std::endl;
                err = current_err;
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
            std::vector<cvm::real> averaged_normal_vector = compute_averaged_border_normal_gradient(
                corresponding_index_in_small_grid);
            for (size_t i = 0; i < nd; i++) {
                potential_value += averaged_normal_vector[i] * relative_position[i] * widths[i];
            }
        }
        data[address(ix)] = potential_value;
    }
};
// In reality, we don't actually need to do that we can just use normal extrapolate, though it's not coherent
void colvargrid_integrate::extrapolate_potential_unweighted() {
    for (std::vector<int> ix = new_index(); index_ok(ix);
         incr(ix)) {
        std::vector<int> corresponding_index_in_small_grid(ix);
        for (size_t i = 0; i < nd; i++) {
            if (!periodic[i]) {
                corresponding_index_in_small_grid[i] = ix[i] - 1;
            }
        }
        std::vector<int> reference_index_in_small_grid(nd, 0);
        bool need_to_extrapolate =
                computation_grid->wrap_to_edge(corresponding_index_in_small_grid, reference_index_in_small_grid);
        std::vector<cvm::real> relative_position(nd, 0);
        cvm::real potential_value = computation_grid->data[computation_grid->address(reference_index_in_small_grid)];
        double only_one = 0;
        std::vector<int> reference_index_in_small_grid_maybe(reference_index_in_small_grid);
        if (need_to_extrapolate) {
            for (size_t i = 0; i < nd; i++) {
                relative_position[i] = corresponding_index_in_small_grid[i] - reference_index_in_small_grid[i];
                only_one += std::abs(relative_position[i]);
                if (relative_position[i] !=0) {
                    reference_index_in_small_grid_maybe[i] -= relative_position[i];
                }
            }
            bool two_or_more = only_one > 1;
            if (!two_or_more){
                potential_value = computation_grid->data[computation_grid->address(reference_index_in_small_grid_maybe)];
                std::vector<cvm::real> averaged_normal_vector = compute_averaged_border_normal_gradient_classic(
                   corresponding_index_in_small_grid);
                for (size_t i = 0; i < nd; i++) {
                    potential_value += averaged_normal_vector[i] * relative_position[i] * 2 * widths[i];
                }
            }
            else {
                potential_value = computation_grid->data[computation_grid->address(reference_index_in_small_grid)];
                std::vector<cvm::real> averaged_normal_vector = compute_averaged_border_normal_gradient(
                   corresponding_index_in_small_grid);
                for (size_t i = 0; i < nd; i++) {
                    potential_value += averaged_normal_vector[i] * relative_position[i] * widths[i];
                }
            }
        }
        data[address(ix)] = potential_value;
        }
};

std::vector<cvm::real> colvargrid_integrate::compute_averaged_border_normal_gradient_classic(const std::vector<int> & ix0) {
    //TODO: maybe do it while computing the div and store it in a vector<vector<real>>>
    std::vector<int> ix(ix0);
    std::vector<cvm::real> averaged_gradient(nd, 0);
    if (nd == 2) {
        std::vector<cvm::real> g00(2,0), g01(2,0), g10(2,0), g11(2,0);
        ix[0] = ix0[0] + 1;
        get_grad(g10, ix);
        ix[1] = ix0[1] + 1;
        get_grad(g11, ix);
        ix[0] = ix0[0];
        get_grad(g01, ix);

        if ((ix0[0] == 0 || ix0[0] == nx[0]) && !periodic[0]) {
            averaged_gradient[0] =  (g10[0] + g00[0] + g11[0] + g01[0]);
        } else if ((ix0[1] == 0 || ix0[1] == nx[1]) && !periodic[1]) {
            averaged_gradient[1] += (g01[1] + g00[1] + g11[1] + g10[1]);
        }
    }
    else if (nd ==3) {
        std::vector<cvm::real> g000(3, 0.0), g001(3, 0.0), g010(3, 0.0), g011(3, 0.0);
        std::vector<cvm::real> g100(3, 0.0), g101(3, 0.0), g110(3, 0.0), g111(3, 0.0);

        // Start at ix0
        ix[0] = ix0[0];
        ix[1] = ix0[1];
        ix[2] = ix0[2];

        // Fill the 8 corners, shifting with +1
        get_grad(g000, ix);            // (0,0,0)

        ix[2] = ix0[2] + 1;
        get_grad(g001, ix);            // (0,0,1)

        ix[1] = ix0[1] + 1;
        ix[2] = ix0[2];
        get_grad(g011, ix);            // (0,1,1)

        ix[2] = ix0[2] - 1 + 1; // or reset properly
        ix[2] = ix0[2];
        get_grad(g010, ix);            // (0,1,0)

        ix[0] = ix0[0] + 1;
        ix[1] = ix0[1];
        ix[2] = ix0[2];
        get_grad(g100, ix);            // (1,0,0)

        ix[2] = ix0[2] + 1;
        get_grad(g101, ix);            // (1,0,1)

        ix[1] = ix0[1] + 1;
        ix[2] = ix0[2];
        get_grad(g111, ix);            // (1,1,1)

        ix[2] = ix0[2] - 1 + 1; // ensure reset
        ix[2] = ix0[2];
        get_grad(g110, ix);            // (1,1,0)

        if ((ix0[0] == 0 || ix0[0] == nx[0]) && !periodic[0]) {
            averaged_gradient[0] =  (g100[0] + g000[0] + g101[0] + g001[0] + g110[0] + g010[0] + g111[0] + g011[0])
                    * 0.125;
        } else if ((ix0[1] == 0 || ix0[1] == nx[1]) && !periodic[1]) {

            averaged_gradient[1]=  (g010[1] + g000[1] + g011[1] + g001[1] + g110[1] + g100[1] + g111[1] + g101[1])
                    * 0.125;
        } else if ((ix0[2] == 0 || ix0[2] == nx[2]) && !periodic[2]) {

            averaged_gradient[2]= (g001[2] - g000[2] + g011[2] - g010[2] + g101[2] - g100[2] + g111[2] - g110[2])
                    * 0.125 ;
        }
    }
    return averaged_gradient;
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

void colvargrid_integrate::extrapolate_data() {
    //TODO: add the count thresholds' calculation ?
    bool converged = false;
    size_t step = 1;
    std::vector<cvm::real> interpolated(gradients->data);
    std::vector<size_t> known(gradients->nt, nt);
    for (std::vector<int> ix = new_index(); gradients->index_ok(ix); incr(ix)) {
        if (gradients->samples->value(ix) > lower_threshold_count) {
            known[gradients->address(ix)] = step;
            interpolated[gradients->address(ix)] = gradients->data[address(ix)];
        }
    }
    std::vector<cvm::real> neighbor_grad(nd, 0);
    cvm::real weight_total = 0.;
    cvm::real weight_neighbor = 0.;
    size_t known_neighbors = 0;
    while (!converged) {
        step++;
        converged = true;
        for (std::vector<int> ix = new_index(); gradients->index_ok(ix); incr(ix)) {
            if (known[gradients->address(ix)] > step) {
                known_neighbors = 0;
                std::vector<cvm::real> interpolated_gradient(nd, 0);
                weight_total = 0;
                for (size_t neighbor = 0; neighbor < laplacian_stencil.size(); neighbor++) {
                    std::vector<int> neighbor_coordinates(ix);
                    for (size_t d = 0; d < nd; d++) {
                        neighbor_coordinates[d] += laplacian_stencil[neighbor][d];
                    }
                    if (known[gradients->address(neighbor_coordinates)] < step) {
                        known_neighbors++;
                        get_grad(neighbor_grad, neighbor_coordinates);
                        weight_neighbor = get_regularized_weight(neighbor_coordinates);
                        weight_total += weight_neighbor;
                        for (size_t d = 0; d < nd; d++) {
                            interpolated_gradient[d] += neighbor_grad[d] * weight_neighbor;
                        }
                    }
                }
                if (known_neighbors > 0) {
                    for (size_t d = 0; d < nd; d++) {
                        interpolated[address(ix) * nd + d] = interpolated_gradient[d] / weight_total;
                    }
                    known[gradients->address(ix)] = step;
                } else
                    converged = false;
            }
        }
    }
    gradients->data = interpolated;
    interpolated.clear();
    known.clear();
}
