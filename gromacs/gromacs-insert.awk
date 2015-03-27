# Insert necessary pieces of code for Colvars into the GROMACS source files.

# Insert our new header.
/#include \"gromacs\/utility\/sysinfo.h\"/ {
    print "#include \"gromacs/pulling/colvars_potential.h\""
}

# Insert the the colvars_potential_wrapper() function.
/static void pull_potential_wrapper/ {
    printf "static void colvars_potential_wrapper(FILE *fplog,t_commrec *cr, t_inputrec*ir, matrix box, rvec x[], rvec f[], tensor vir_force, t_mdatoms *mdatoms,gmx_enerdata_t *enerd, gmx_int64_t step, gmx_wallcycle_t wcycle) {\nt_pbc pbc;\nwallcycle_start(wcycle, ewcPULLPOT);\nset_pbc(&pbc, ir->ePBC, box);\nenerd->term[F_COM_PULL] += colvars_potential(ir, mdatoms, &pbc, step, x, f, vir_force);\nwallcycle_stop(wcycle, ewcPULLPOT);\n}\n\n"
}

# Insert the calls to colvars_potential_wrapper().
/inputrec->bPull && inputrec->pull->bPotential/ {
    printf "/* Colvars operation is currently defined by userint1*/\nif (0 != inputrec->userint1) colvars_potential_wrapper(fplog, cr, inputrec, box, x, f, vir_force, mdatoms, enerd, step, wcycle);\n\n" 
}

# Keep the already existing code.
{
    print $0
}
