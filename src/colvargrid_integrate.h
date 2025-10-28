#ifndef COLVARGRID_INTEGRATE_H
#define COLVARGRID_INTEGRATE_H

#include <iostream>

#include "colvargrid.h"

/// Integrate (1D, 2D or 3D) gradients
class colvargrid_integrate : public colvar_grid_scalar {
public:
    colvargrid_integrate() {
    };

    virtual ~colvargrid_integrate() {
    }


    /// Constructor from a vector of colvars + gradient grid
    colvargrid_integrate(std::vector<colvar *> &colvars,
                         std::shared_ptr<colvar_grid_gradient> gradients);

    /// Constructor from a gradient grid (for processing grid files without a Colvars config)
    colvargrid_integrate(std::shared_ptr<colvar_grid_gradient> gradients, bool is_weighted = false);

    /// \brief Calculate potential from divergence (in 2D); return number of steps
    int integrate(const int itmax, const cvm::real &tol, cvm::real &err, bool verbose = true);

    /// \brief Update matrix containing divergence and boundary conditions
    /// based on new gradient point value, in neighboring bins
    void update_div_neighbors(const std::vector<int> &ix);

    /// \brief Update matrix containing divergence and boundary conditions
    /// called by update_div_neighbors and by colvarbias_abf::adiabatic_reweighting_update_gradient_pmf
    void update_weighted_div_local(const std::vector<int> &ix);

    void update_div_local(const std::vector<int> &ix0);

    /// \brief Set matrix containing divergence and boundary conditions
    /// based on complete gradient grid
    void set_div();

    void set_weighted_div();


    /// \brief Add constant to potential so that its minimum value is zero
    /// Useful e.g. for output
    inline void set_zero_minimum() {
        add_constant(-1.0 * minimum_value());
    }

    /// \brief Flag requesting the use of a smoothed version of the gradient (default: false)
    bool b_smoothed;

    // \brief For testing purposes only: print the different stencils computed in prepare_laplacian_calculation.
    void print_laplacian_preparations();

    // \brief Computes all the relative positions of objects necessary to calculate the laplacian at a specific point
    void prepare_laplacian_necessary_stencils();

    // \brief Computes all the relative positions to calculate the divergence at a specific point
    void prepare_divergence_necessary_stencils();

    void prepare_calculations();

    /// Array holding divergence + boundary terms (modified Neumann) if not periodic
    std::vector<cvm::real> laplacian_matrix_test;

protected:
    bool weighted = false;
    bool precompute = true;

    colvar_grid_scalar *computation_grid = new colvar_grid_scalar();
    std::vector<cvm::real> div_border_supplement;

    std::vector<int> computation_nx;
    std::vector<int> computation_nxc;

    size_t computation_nt;
    bool need_to_extrapolate_solution = false;
    // Reference to gradient grid
    std::shared_ptr<colvar_grid_gradient> gradients;
    std::vector<cvm::real> divergence;


    // Scalar grid containing interpolated weights, same mesh as FES and Laplacian
    // Stored as a flat vector like the divergence
    std::vector<cvm::real> weights;
    std::vector<cvm::real> regularized_weights;
    std::vector<cvm::real> laplacian_coefficients;
    std::vector<size_t> sorted_counts;

    size_t sum_count;
    // max and min count to regularize F
    size_t max_count_F = 1;
    size_t min_count_F = 0;
    // max and min count to regularize the weights
    float lambda_max = 0.3;
    float lambda_min = 0.3;
    size_t upper_threshold_count = 1;
    size_t lower_threshold_count = 1;
    size_t m_num_threads = 1;

    // Get G at a specific point where G is the gradient F if there is enough observation else it's F multiplied by a coefficient < 1
    void get_regularized_F(std::vector<cvm::real> &F, std::vector<int> &ix);

    // Get weight regularized by a lower and upper threshold and a ramp in between
    cvm::real get_regularized_weight(std::vector<int> &ix);

    // positions of the points in the stencil relative to the stencil center
    std::vector<std::vector<int> > laplacian_stencil;
    // positions of the weights relative to each stencil point to take into account in the weighted laplacian
    std::vector<std::vector<std::vector<int> > > weight_stencil;
    // Coefficient of each point in the stencil
    std::vector<cvm::real> weight_counts;
    // for each point in the stencil tells if it is also included in the classical laplacian stencil and what is its
    // coefficient
    std::vector<cvm::real> neighbor_in_classic_laplacian_stencil;
    // relative coordinates (in the data grid) of the points to take into account in the divergence calculation
    std::vector<std::vector<int> > surrounding_points_relative_positions;
    //   std::vector<cvm::real> inv_lap_diag; // Inverse of the diagonal of the Laplacian; for conditioning

    /// Obtain the gradient vector at given location ix, if available
    /// ix needs to be wrapped beforehand in PBC, if necessary
    /// Returns the sample count in given bin if available.
    size_t get_grad(std::vector<cvm::real> &g, std::vector<int> ix);

    /// \brief Solve linear system based on CG, valid for symmetric matrices only
    /// atimes : left multiplication by LHS symmetric matrix
    /// b : RHS of equation
    /// x : initial guess for the solution; output is solution
    /// tol : convergence criterion
    /// itmax : max number of iterations
    /// iter (out) : final iteration
    /// err (out) : final error value
    void nr_linbcg_sym(const bool weighted, const std::vector<cvm::real> &b, std::vector<cvm::real> &x,
                       const cvm::real &tol, const int itmax, int &iter, cvm::real &err);

    void optimize_adam(bool weighted, const std::vector<cvm::real> &b, std::vector<cvm::real> &x,
                       const cvm::real &tol, const int itmax, int &iter, cvm::real &err);

    /// l2 norm of a vector
    cvm::real l2norm(const std::vector<cvm::real> &x);

    /// @brief Multiplication by sparse matrix representing Lagrangian (or its transpose)
    /// @param x scalar field, discretized on grid, flattened
    /// @param r (out) discrete Laplacian of x
    void laplacian(const std::vector<cvm::real> &x, std::vector<cvm::real> &r);

    template<bool initialize_div_supplement>
    void laplacian_weighted(const std::vector<cvm::real> &x, std::vector<cvm::real> &r);

    /// Computes the line result of the weighted laplacian matrix multiplied by x and stores it in the appropriate line
    /// in r. Uses the coefficients that must be precomputed.
    template<bool initialize_div_supplement>
    void linewise_laplacian_weighted_precomputed(const std::vector<cvm::real> &x, std::vector<cvm::real> &r,
                                                 size_t grid_address);
    /// Computes the line result of the weighted laplacian matrix multiplied by x and stores it in the appropriate line
    /// in r. Calculates laplacian coefficients on the fly.
    template<bool initialize_div_supplement>
    void linewise_laplacian_weighted_otf(const std::vector<cvm::real> &x, std::vector<cvm::real> &r,
                                         size_t grid_address);
    typedef void (colvargrid_integrate::*func_pointer)(const std::vector<cvm::real> &, std::vector<cvm::real> &, size_t);
    func_pointer linewise_laplacian_weighted;
    /// Compute gradient of whole potential grid by finite difference
    // void compute_grad(const std::vector<cvm::real> &A, std::vector<cvm::real> &G);


    //   /// Inversion of preconditioner matrix
    //   void asolve(const std::vector<cvm::real> &b, std::vector<cvm::real> &x);

    /// Converts a number in base 3.
    std::string convert_base_three(int n);

    std::string convert_base_two(int n, size_t length);

    /// Computes the normal border gradient for ghost points calculations in the weighted scheme.
    std::vector<cvm::real> compute_averaged_border_normal_gradient(std::vector<int> virtual_point_coordinates);

    /// Computes the normal border gradient for ghost points calculations in the unweighted scheme, i.e. like Long Chen.
    std::vector<cvm::real> compute_averaged_border_normal_gradient_classic(const std::vector<int> &ix0);

    template<typename T>
    typename std::vector<T>::iterator insert_into_sorted_list(std::vector<T> &sortedList, const T &value);

    inline void reverse(std::string::iterator, std::string::iterator);

    /// From the smaller resolution grid, use order 1taylor expansion to find the value on a grid with 2 more cell in each
    /// dimension. This means, we assigned the ghost points values which are the ones we assigned during the laplacian
    /// calculation.
    void extrapolate_potential();

    // Same as extrapolate potential, but we assign the values we assumed in the unweighted laplacian scheme. For corner
    // points also use order 1 Taylor expansion.
    void extrapolate_potential_unweighted();

    /// \brief From the cells with estimate of the gradient of the free energy propagates the information through
    /// interpolation to other cells where data is lacking
    void extrapolate_data();

    /// \brief Initialize computation_nx based on nx and periodic boundaries
    inline int init_computation_nx_nt() {
        computation_nx.resize(nd);
        computation_nt = 1;
        computation_nxc.resize(nd);
        for (size_t i = 0; i < nd; i++) {
            if (periodic[i]) {
                computation_nx[i] = nx[i];
            } else {
                computation_nx[i] = nx[i] - 1; // One less point for non-periodic dimensions
            }
            computation_nt *= computation_nx[i];
            computation_nxc[i] = computation_nt;
        }
#ifdef _OPENMP
        m_num_threads = cvm::proxy->smp_num_threads();
#else
        if (m_num_threads > 1) {
            return cvm::error(
                "Multi-threading requested in weighted integrator, which is not supported by this build.\n");
        }
#endif
        std::cout << m_num_threads << " : nb threads" << std::endl;
        return 0;
    }
};
#endif
