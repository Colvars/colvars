#include <iostream>
#include <fstream>
#include <sys/stat.h>

#include "colvargrid.h"
#include "colvargrid_integrate.h"
#include "colvarproxy.h"

void saveVectorToCSV(const std::vector<cvm::real> &vec, const std::string &filename)
{
    std::ofstream file(filename);
    if (!file)
    {
        std::cerr << "Error opening file\n";
        return;
    }

    for (size_t i = 0; i < vec.size(); ++i)
    {
        file << vec[i];
        if (i != vec.size() - 1)
            file << ","; // Separate values with commas
    }
    file.close();
}

int main(int argc, char *argv[])
{
    bool weighted = true;
    bool save_divergence = false; // For testing: need to uncomment lines and put divergence back in public
    int itmax = 1000;
    cvm::real err;
    cvm::real tol = 2e-3;
    if (argc < 2)
    {
        std::cerr << "\n\nOne argument needed: gradient multicol file name.\n";
        return 1;
    }
    colvarproxy *proxy = new colvarproxy();
    proxy->colvars = new colvarmodule(proxy); // This could be omitted if we used the colvarproxy_stub class

    std::string gradfile(argv[1]);
    std::string countfile;
    std::shared_ptr<colvar_grid_count> count_ptr;

    // Look for matching count file
    if (argc == 2)
    {
        size_t pos = gradfile.rfind(std::string(".czar.grad"));
        if (pos != std::string::npos)
        {
            countfile = gradfile.substr(0, pos) + ".zcount";
        }
        else
        {
            pos = gradfile.rfind(std::string(".grad"));
            if (pos != std::string::npos)
            {
                countfile = gradfile.substr(0, pos) + ".count";
            }
        }
    }
    else if (argc > 2)
    {
        countfile = argv[2];
    }

    if (countfile.size())
    {
        struct stat buffer;
        if (stat(countfile.c_str(), &buffer) == 0)
        {
            std::cout << "Found associated count file " << countfile << ", reading...\n";
            count_ptr.reset(new colvar_grid_count(countfile));
            if (!count_ptr || count_ptr->nd == 0)
            {
                // catch constructor failure
                cvm::error("Error reading count grid.");
                return cvm::get_error();
            }
        }
    }

    count_ptr->write_multicol("counts.dat");

    std::cout << "Reading gradient file " << gradfile << std::endl;
    std::shared_ptr<colvar_grid_gradient> grad_ptr = std::make_shared<colvar_grid_gradient>(gradfile, count_ptr);
    // std::shared_ptr<colvar_grid_gradient> grad_ptr = std::make_shared<colvar_grid_gradient>(gradfile);

    if (!grad_ptr || grad_ptr->nd == 0)
    {
        // catch constructor failure
        cvm::error("Error reading gradient grid.");
        return cvm::get_error();
    }

    grad_ptr->write_multicol("gradient_in.dat");

    colvargrid_integrate fes(grad_ptr, weighted);

    // ####################### TESTING ##################################
    // fes.prepare_laplacian_necessary_stencils();
    // fes.print_laplacian_preparations();
    // If want to test by saving div then put div_supplement and divergence in public.
    // if (save_divergence) {
    //     if (!weighted) {
    //         fes.set_div();
    //     } else {
    //         fes.prepare_calculations();
    //         fes.set_weighted_div();
    //         // fes.laplacian_weighted<true>(fes.divergence, fes.data);
    //         for (size_t i = 0; i < itmax; i++) {
    //             fes.divergence[i] = fes.divergence[i];
    //             // if (fes.div_border_supplement[i] > tol) {
    //                 // std::cout << "ola ";
    //                 // std::cout << fes.div_border_supplement[i] << std::endl;
    //             }
    //         }
    //     }
    // colvar_grid_scalar div(fes);
    // div.data = fes.divergence;
    // div.nx = fes.computation_nx;
    // std::cout << div.nx[0] << " " << div.nx[1] << std::endl;
    // div.nt = fes.computation_nt;
    // div.nxc = fes.computation_nxc;
    // div.write_multicol("divergence.dat");
    // std::cout << "\nWriting divergence in multicol format to divergence.dat";
    // saveVectorToCSV(fes.laplacian_matrix_test, "laplacian.csv");
    // }
    // std::vector<cvm::real> laplacian_matrix (fes.computation_grid->nt, 0);
    // std::vector<cvm::real> test_vector (fes.computation_grid->nt, 1);
    // // std::vector<cvm::real> complete_div (fes.computation_grid->nt, 0);
    // // //
    // fes.laplacian_weighted<true>(test_vector, laplacian_matrix);
    // for (int i = 0; i < fes.computation_grid->nt; i++){
    // complete_div[i] = fes.divergence[i] + fes.div_border_supplement[i];
    // }
    // std::cout << fes.divergence.size() << " " << fes.div_border_supplement.size() << std::endl;
    // saveVectorToCSV(fes.divergence, "divergence.csv");

    fes.integrate(itmax, tol, err, true);
    fes.set_zero_minimum();
    if (fes.num_variables() < 3)
    {
        if (weighted)
        {
            fes.write_multicol(std::string(gradfile + ".int"), "integrated fes");
            std::cout << "\nWriting integrated fes in multicol format to " + gradfile + ".int\n";
        }
        else
        {
            fes.write_multicol(std::string(gradfile + ".int"), "integrated fes");
            std::cout << "\nWriting integrated fes in multicol format to " + gradfile + ".int\n";
        }
    }
    else
    {
        // Write 3D grids to more convenient DX format
        std::cout << "\nWriting integrated free energy in OpenDX format to " + gradfile + ".int.dx\n";
        fes.write_opendx(std::string(gradfile + ".int.dx"), "integrated free energy");
    }
    saveVectorToCSV(fes.laplacian_matrix_test, "laplacian.csv");

    delete proxy;
    return 0;
}
