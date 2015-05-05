// Jeff Comer's tests to see if he can link GROMACS and Colvars
#ifndef GMX_COLVARS_COLVARS_POTENTIAL
#define GMX_COLVARS_COLVARS_POTENTIAL

#ifdef __cplusplus
extern "C" {
#endif
real colvars_potential(t_inputrec *ir, t_mdatoms *md, t_pbc *pbc,
			 gmx_int64_t step, rvec *x, rvec *f, tensor vir);
#ifdef __cplusplus
}
#endif

// Implemented in 'colvarproxy_gromacs.cpp'
#endif
