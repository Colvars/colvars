#include <iostream>
#include <fstream>
#include <sys/stat.h>

#include "colvargrid.h"
#include "colvargrid_integrate.h"
#include "colvarproxy.h"


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

    std::cout << "Reading gradient file " << gradfile << std::endl;
    std::shared_ptr<colvar_grid_gradient> grad_ptr = std::make_shared<colvar_grid_gradient>(gradfile, count_ptr);

    if (!grad_ptr || grad_ptr->nd == 0)
    {
        // catch constructor failure
        cvm::error("Error reading gradient grid.");
        return cvm::get_error();
    }

    grad_ptr->write_multicol("gradient_in.dat");

    colvargrid_integrate fes(grad_ptr, weighted);

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

    delete proxy;
    return 0;
}
