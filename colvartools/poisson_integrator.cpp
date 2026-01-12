#include <iostream>
#include <fstream>

#include "colvargrid.h"
#include "colvargrid_integrate.h"
#include "colvarproxy.h"
#include "CLI11.hpp"

int main(int argc, char *argv[])
{
    bool weighted = false;
    bool debug = false;
    int itmax = 1000;
    cvm::real err;
    cvm::real tol = 2e-3;

    std::string gradfile;
    std::string countfile;

    CLI::App app("Colvars gradient grid integrator");

    app.add_option("gradfile", gradfile, "Gradient multicol file name")
        ->required()
        ->check(CLI::ExistingFile);

    app.add_option("countfile", countfile, "Count file name");
    app.add_flag("--weighted,!--unweighted", weighted, "Enable or disable weighting")
        ->default_val(weighted);

    app.add_option("--tol", tol, "Tolerance value");
    app.add_option("--itmax", itmax, "Maximum iterations");
    app.add_flag("--debug", debug, "Write debugging output and files");

    CLI11_PARSE(app, argc, argv);

    // (Moved after parsing so we don't allocate memory if help/error is printed)
    colvarproxy *proxy = new colvarproxy();
    proxy->colvars = new colvarmodule(proxy);
    std::shared_ptr<colvar_grid_count> count_ptr;

    // Deduce countfile if not provided by user
    if (countfile.empty())
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
            if (debug) count_ptr->write_multicol("counts.dat");
        }
    }

    std::cout << "Reading gradient file " << gradfile << std::endl;
    std::shared_ptr<colvar_grid_gradient> grad_ptr = std::make_shared<colvar_grid_gradient>(gradfile, count_ptr);

    if (!grad_ptr || grad_ptr->nd == 0 || grad_ptr->nd != grad_ptr->multiplicity() || cvm::get_error() != COLVARS_OK)
    {
        // catch constructor failure
        cvm::error("Error reading gradient grid.");
        return cvm::get_error();
    }

    if (debug) grad_ptr->write_multicol("gradient_in.dat");

    colvargrid_integrate fes(grad_ptr, weighted);
    fes.integrate(itmax, tol, err, true);
    fes.set_zero_minimum();
    if (fes.num_variables() < 3)
    {
        std::cout << "\nWriting integrated fes in multicol format to " + gradfile + ".int\n";
        fes.write_multicol(std::string(gradfile + ".int"), "integrated fes");
    } else {
        // Write 3D grids to more convenient DX format
        std::cout << "\nWriting integrated free energy in OpenDX format to " + gradfile + ".int.dx\n";
        fes.write_opendx(std::string(gradfile + ".int.dx"), "integrated free energy");
    }

    delete proxy;
    return 0;
}
