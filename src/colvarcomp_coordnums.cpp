// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include "colvarmodule.h"
#include "colvaratoms.h"
#include "colvarvalue.h"
#include "colvar.h"
#include "colvarcomp.h"
#include "colvarcomp_coordnums.h"


colvar::coordnum::coordnum()
{
  set_function_type("coordNum");
  x.type(colvarvalue::type_scalar);
  cvm::real const r0 = cvm::main()->proxy->angstrom_to_internal(4.0);
  r0_vec = {r0, r0, r0};
  // Default upper boundary not yet known
}


int colvar::coordnum::init(std::string const &conf)
{
  int error_code = cvc::init(conf);

  group1 = parse_group(conf, "group1");

  if (!group1) {
    return error_code | COLVARS_INPUT_ERROR;
  }

  if (group1->b_dummy) {
    error_code |= cvm::error("Error: group1 may not be a dummy atom\n", COLVARS_INPUT_ERROR);
  }

  if (function_type() != "selfCoordNum") {

    group2 = parse_group(conf, "group2");
    if (!group2) {
      return error_code | COLVARS_INPUT_ERROR;
    }

    if (int atom_number = cvm::atom_group::overlap(*group1, *group2)) {
      error_code |= cvm::error("Error: group1 and group2 share a common atom (number: " +
                                   cvm::to_str(atom_number) + ")\n",
                               COLVARS_INPUT_ERROR);
    }

    if (function_type() == "coordNum") {
      get_keyval(conf, "group1CenterOnly", b_group1_center_only, group2->b_dummy);
      get_keyval(conf, "group2CenterOnly", b_group2_center_only, group2->b_dummy);
    }

    if (function_type() == "groupCoord") {
      // In groupCoord, these flags are hard-coded
      b_group1_center_only = true;
      b_group2_center_only = true;
    }

    size_t const group1_num_coords = b_group1_center_only ? 1 : group1->size();
    size_t const group2_num_coords = b_group2_center_only ? 1 : group2->size();

    num_pairs = group1_num_coords * group2_num_coords;

  } else {

    // selfCoordNum case
    num_pairs = (group1->size() * (group1->size() - 1)) / 2;
  }

  init_scalar_boundaries(0.0, num_pairs);

  // Get the default value from r0_vec to report it
  cvm::real r0 = r0_vec[0];
  bool const b_redefined_cutoff = get_keyval(conf, "cutoff", r0, r0);

  if (get_keyval(conf, "cutoff3", r0_vec, r0_vec)) {
    if (b_redefined_cutoff) {
      error_code |=
          cvm::error("Error: cannot specify \"cutoff\" and \"cutoff3\" at the same time.\n",
                     COLVARS_INPUT_ERROR);
    }

    // remove meaningless negative signs
    if (r0_vec.x < 0.0) r0_vec.x *= -1.0;
    if (r0_vec.y < 0.0) r0_vec.y *= -1.0;
    if (r0_vec.z < 0.0) r0_vec.z *= -1.0;
  } else {
    if (b_redefined_cutoff) {
      r0_vec = {r0, r0, r0};
    }
  }

  get_keyval(conf, "expNumer", en, en);
  get_keyval(conf, "expDenom", ed, ed);

  if ( (en%2) || (ed%2) ) {
    error_code |= cvm::error("Error: odd exponent(s) provided, can only use even ones.\n",
                             COLVARS_INPUT_ERROR);
  }

  if ( (en <= 0) || (ed <= 0) ) {
    error_code |= cvm::error("Error: negative exponent(s) provided.\n",
                             COLVARS_INPUT_ERROR);
  }

  if (!is_enabled(f_cvc_pbc_minimum_image)) {
    cvm::log("Warning: only minimum-image distances are used by this variable.\n");
  }

  if (function_type() != "groupCoord") {
    // All coordNum variables may benefit from a pairlist, except groupCoord
    get_keyval(conf, "tolerance", tolerance, tolerance);
    if (tolerance > 0) {
      cvm::main()->cite_feature("coordNum pairlist");
      get_keyval(conf, "pairListFrequency", pairlist_freq, pairlist_freq);
      if ( ! (pairlist_freq > 0) ) {
        return cvm::error("Error: non-positive pairlistfrequency provided.\n",
                          COLVARS_INPUT_ERROR);
        // return and do not allocate the pairlists below
      }
      pairlist = new bool[num_pairs];
      for (size_t i = 0; i < num_pairs; i++) pairlist[i] = true;
    }
  }

  return error_code;
}


colvar::coordnum::~coordnum()
{
  if (pairlist) {
    delete [] pairlist;
  }
}


template<int flags> void colvar::coordnum::main_loop(bool **pairlist_elem)
{
  const cvm::rvector inv_r0_vec{
    1.0 / r0_vec.x,
    1.0 / r0_vec.y,
    1.0 / r0_vec.z};
  cvm::rvector const inv_r0sq_vec{
    inv_r0_vec.x*inv_r0_vec.x,
    inv_r0_vec.y*inv_r0_vec.y,
    inv_r0_vec.z*inv_r0_vec.z};
  if (b_group2_center_only) {
    const cvm::atom_pos group2_com = group2->center_of_mass();
    cvm::rvector group2_com_grad(0, 0, 0);
    for (size_t i = 0; i < group1->size(); ++i) {
      x.real_value += switching_function<flags>(inv_r0_vec,
                                                inv_r0sq_vec, en, ed,
                                                group1->pos_x(i),
                                                group1->pos_y(i),
                                                group1->pos_z(i),
                                                group2_com.x,
                                                group2_com.y,
                                                group2_com.z,
                                                group1->grad_x(i),
                                                group1->grad_y(i),
                                                group1->grad_z(i),
                                                group2_com_grad.x,
                                                group2_com_grad.y,
                                                group2_com_grad.z,
                                                pairlist_elem,
                                                tolerance);
    }
    group2->set_weighted_gradient(group2_com_grad);
  } else {
    for (size_t i = 0; i < group1->size(); ++i) {
      for (size_t j = 0; j < group2->size(); ++j) {
        x.real_value += switching_function<flags>(inv_r0_vec,
                                                  inv_r0sq_vec, en, ed,
                                                  group1->pos_x(i),
                                                  group1->pos_y(i),
                                                  group1->pos_z(i),
                                                  group2->pos_x(j),
                                                  group2->pos_y(j),
                                                  group2->pos_z(j),
                                                  group1->grad_x(i),
                                                  group1->grad_y(i),
                                                  group1->grad_z(i),
                                                  group2->grad_x(j),
                                                  group2->grad_y(j),
                                                  group2->grad_z(j),
                                                  pairlist_elem,
                                                  tolerance);
      }
    }
  }
}


template<int compute_flags> int colvar::coordnum::compute_coordnum()
{
  bool const use_pairlist = (pairlist != NULL);
  bool const rebuild_pairlist = (pairlist != NULL) &&
    (cvm::step_relative() % pairlist_freq == 0);

  bool *pairlist_elem = use_pairlist ? pairlist : NULL;

  if (use_pairlist) {
    if (rebuild_pairlist) {
      int const flags = compute_flags | ef_use_pairlist | ef_rebuild_pairlist;
      main_loop<flags>(&pairlist_elem);
    } else {
      int const flags = compute_flags | ef_use_pairlist;
      main_loop<flags>(&pairlist_elem);
    }
  } else {
    int const flags = compute_flags;
    main_loop<flags>(NULL);
  }

  return COLVARS_OK;
}


void colvar::coordnum::calc_value()
{
  x.real_value = 0.0;
  if (is_enabled(f_cvc_gradient)) {
    compute_coordnum<ef_gradients>();
  } else {
    compute_coordnum<ef_null>();
  }
}


void colvar::coordnum::calc_gradients()
{
  // Gradients are computed by calc_value() if f_cvc_gradients is enabled
}



// h_bond member functions

colvar::h_bond::h_bond()
{
  cvm::real const r0 = cvm::main()->proxy->angstrom_to_internal(3.3);
  r0_vec = {r0, r0, r0};
  set_function_type("hBond");
  x.type(colvarvalue::type_scalar);
  init_scalar_boundaries(0.0, 1.0);
}


int colvar::h_bond::init(std::string const &conf)
{
  int error_code = cvc::init(conf);

  if (cvm::debug())
    cvm::log("Initializing h_bond object.\n");

  set_function_type("hBond");
  x.type(colvarvalue::type_scalar);
  init_scalar_boundaries(0.0, 1.0);

  int a_num = -1, d_num = -1;
  get_keyval(conf, "acceptor", a_num, a_num);
  get_keyval(conf, "donor",    d_num, a_num);

  if ( (a_num == -1) || (d_num == -1) ) {
    error_code |= cvm::error("Error: either acceptor or donor undefined.\n", COLVARS_INPUT_ERROR);
  }

  register_atom_group(new cvm::atom_group);
  {
    colvarproxy* const p = cvm::main()->proxy;
    auto modify_atom = atom_groups[0]->get_atom_modifier();
    modify_atom.add_atom(cvm::atom_group::init_atom_from_proxy(p, a_num));
    modify_atom.add_atom(cvm::atom_group::init_atom_from_proxy(p, d_num));
  }

  cvm::real r0 = r0_vec[0];
  bool const b_redefined_cutoff = get_keyval(conf, "cutoff", r0, r0);
  if (b_redefined_cutoff) {
    r0_vec = {r0, r0, r0};
  }
  get_keyval(conf, "expNumer", en, en);
  get_keyval(conf, "expDenom", ed, ed);

  if ((en % 2) || (ed % 2)) {
    error_code |= cvm::error("Error: odd exponent(s) provided, can only use even ones.\n",
                             COLVARS_INPUT_ERROR);
  }

  if ((en <= 0) || (ed <= 0)) {
    error_code |= cvm::error("Error: negative exponent(s) provided.\n", COLVARS_INPUT_ERROR);
  }

  if (cvm::debug())
    cvm::log("Done initializing h_bond object.\n");

  return error_code;
}

colvar::h_bond::h_bond(cvm::atom_group::simple_atom const &acceptor,
                       cvm::atom_group::simple_atom const &donor,
                       cvm::real r0_i, int en_i, int ed_i)
  : h_bond()
{
  r0_vec = {r0_i, r0_i, r0_i};
  en = en_i;
  ed = ed_i;
  register_atom_group(new cvm::atom_group);
  auto modify_atom = atom_groups[0]->get_atom_modifier();
  modify_atom.add_atom(acceptor);
  modify_atom.add_atom(donor);
}


void colvar::h_bond::calc_value()
{
  int const flags = coordnum::ef_null;
  cvm::rvector G1, G2;
  const cvm::atom_pos A1{atom_groups[0]->pos_x(0),
                         atom_groups[0]->pos_y(0),
                         atom_groups[0]->pos_z(0)};
  const cvm::atom_pos A2{atom_groups[0]->pos_x(1),
                         atom_groups[0]->pos_y(1),
                         atom_groups[0]->pos_z(1)};

  const cvm::rvector inv_r0_vec{
    1.0 / r0_vec.x,
    1.0 / r0_vec.y,
    1.0 / r0_vec.z
  };
  cvm::rvector const inv_r0sq_vec{
    inv_r0_vec.x * inv_r0_vec.x,
    inv_r0_vec.y * inv_r0_vec.y,
    inv_r0_vec.z * inv_r0_vec.z
  };

  x.real_value =
    coordnum::switching_function<flags>(inv_r0_vec, inv_r0sq_vec, en, ed,
                                        atom_groups[0]->pos_x(0),
                                        atom_groups[0]->pos_y(0),
                                        atom_groups[0]->pos_z(0),
                                        atom_groups[0]->pos_x(1),
                                        atom_groups[0]->pos_y(1),
                                        atom_groups[0]->pos_z(1),
                                        atom_groups[0]->grad_x(0),
                                        atom_groups[0]->grad_y(0),
                                        atom_groups[0]->grad_z(0),
                                        atom_groups[0]->grad_x(1),
                                        atom_groups[0]->grad_y(1),
                                        atom_groups[0]->grad_z(1),
                                        NULL, 0.0);
  // Skip the gradient
}


void colvar::h_bond::calc_gradients()
{
  int const flags = coordnum::ef_gradients;
  const cvm::rvector inv_r0_vec{
    1.0 / r0_vec.x,
    1.0 / r0_vec.y,
    1.0 / r0_vec.z
  };
  cvm::rvector const inv_r0sq_vec{
    inv_r0_vec.x*inv_r0_vec.x,
    inv_r0_vec.y*inv_r0_vec.y,
    inv_r0_vec.z*inv_r0_vec.z
  };
  coordnum::switching_function<flags>(inv_r0_vec, inv_r0sq_vec, en, ed,
                                      atom_groups[0]->pos_x(0),
                                      atom_groups[0]->pos_y(0),
                                      atom_groups[0]->pos_z(0),
                                      atom_groups[0]->pos_x(1),
                                      atom_groups[0]->pos_y(1),
                                      atom_groups[0]->pos_z(1),
                                      atom_groups[0]->grad_x(0),
                                      atom_groups[0]->grad_y(0),
                                      atom_groups[0]->grad_z(0),
                                      atom_groups[0]->grad_x(1),
                                      atom_groups[0]->grad_y(1),
                                      atom_groups[0]->grad_z(1),
                                      NULL, 0.0);
}


colvar::selfcoordnum::selfcoordnum()
{
  set_function_type("selfCoordNum");
}


template<int compute_flags> int colvar::selfcoordnum::compute_selfcoordnum()
{
  bool const use_pairlist = (pairlist != NULL);
  bool const rebuild_pairlist = (pairlist != NULL) &&
    (cvm::step_relative() % pairlist_freq == 0);

  bool *pairlist_elem = use_pairlist ? pairlist : NULL;
  size_t i = 0, j = 0;
  size_t const n = group1->size();

  // Always isotropic (TODO: enable the ellipsoid?)
#define CALL_KERNEL(flags) do {                         \
  const cvm::rvector inv_r0_vec{                        \
    1.0 / r0_vec.x,                                     \
    1.0 / r0_vec.y,                                     \
    1.0 / r0_vec.z                                      \
  };                                                    \
  cvm::rvector const inv_r0sq_vec{                      \
    inv_r0_vec.x*inv_r0_vec.x,                          \
    inv_r0_vec.y*inv_r0_vec.y,                          \
    inv_r0_vec.z*inv_r0_vec.z                           \
  };                                                    \
  for (i = 0; i < n - 1; i++) {                         \
    for (j = i + 1; j < n; j++) {                       \
      x.real_value +=                                   \
        coordnum::switching_function<flags>(            \
          inv_r0_vec, inv_r0sq_vec, en, ed,             \
          group1->pos_x(i),                             \
          group1->pos_y(i),                             \
          group1->pos_z(i),                             \
          group1->pos_x(j),                             \
          group1->pos_y(j),                             \
          group1->pos_z(j),                             \
          group1->grad_x(i),                            \
          group1->grad_y(i),                            \
          group1->grad_z(i),                            \
          group1->grad_x(j),                            \
          group1->grad_y(j),                            \
          group1->grad_z(j),                            \
          &pairlist_elem, tolerance);                   \
    }                                                   \
  }                                                     \
} while (0);

  if (use_pairlist) {
    if (rebuild_pairlist) {
      int const flags = compute_flags | coordnum::ef_use_pairlist |
        coordnum::ef_rebuild_pairlist;
      CALL_KERNEL(flags);
    } else {
      int const flags = compute_flags | coordnum::ef_use_pairlist;
      CALL_KERNEL(flags);
    }
  } else { // if (use_pairlist) {
    int const flags = compute_flags | coordnum::ef_null;
      CALL_KERNEL(flags);
  }
#undef CALL_KERNEL
  return COLVARS_OK;
}


void colvar::selfcoordnum::calc_value()
{
  x.real_value = 0.0;
  if (is_enabled(f_cvc_gradient)) {
    compute_selfcoordnum<coordnum::ef_gradients>();
  } else {
    compute_selfcoordnum<coordnum::ef_null>();
  }
}


void colvar::selfcoordnum::calc_gradients()
{
  // Gradients are computed by calc_value() if f_cvc_gradients is enabled
}



colvar::groupcoordnum::groupcoordnum()
{
  set_function_type("groupCoord");
}


void colvar::groupcoordnum::calc_value()
{
  const cvm::atom_pos A1 = group1->center_of_mass();
  const cvm::atom_pos A2 = group2->center_of_mass();
#define CALL_KERNEL(flags) do {                    \
  const cvm::rvector inv_r0_vec{                   \
    1.0 / r0_vec.x,                                \
    1.0 / r0_vec.y,                                \
    1.0 / r0_vec.z                                 \
  };                                               \
  cvm::rvector const inv_r0sq_vec{                 \
    inv_r0_vec.x*inv_r0_vec.x,                     \
    inv_r0_vec.y*inv_r0_vec.y,                     \
    inv_r0_vec.z*inv_r0_vec.z                      \
  };                                               \
  cvm::rvector G1, G2;                             \
  x.real_value = coordnum::switching_function<flags>(inv_r0_vec, inv_r0sq_vec, en, ed, \
                                                     A1.x, A1.y, A1.z, \
                                                     A2.x, A2.y, A2.z, \
                                                     G1.x, G1.y, G1.z, \
                                                     G2.x, G2.y, G2.z, NULL, 0.0); \
} while (0);
  CALL_KERNEL(coordnum::ef_null);
#undef CALL_KERNEL
}


void colvar::groupcoordnum::calc_gradients()
{
  const cvm::atom_pos A1 = group1->center_of_mass();
  const cvm::atom_pos A2 = group2->center_of_mass();
  cvm::rvector G1(0, 0, 0), G2(0, 0, 0);
#define CALL_KERNEL(flags) do { \
  const cvm::rvector inv_r0_vec{              \
    1.0 / r0_vec.x,                           \
    1.0 / r0_vec.y,                           \
    1.0 / r0_vec.z                            \
  };                                          \
  cvm::rvector const inv_r0sq_vec{            \
    inv_r0_vec.x*inv_r0_vec.x,                \
    inv_r0_vec.y*inv_r0_vec.y,                \
    inv_r0_vec.z*inv_r0_vec.z                 \
  };                                          \
  coordnum::switching_function<flags>(inv_r0_vec, inv_r0sq_vec, en, ed, \
                                      A1.x, A1.y, A1.z, \
                                      A2.x, A2.y, A2.z, \
                                      G1.x, G1.y, G1.z, \
                                      G2.x, G2.y, G2.z, NULL, 0.0); \
} while (0);
  int const flags = coordnum::ef_gradients;
  CALL_KERNEL(flags);
  group1->set_weighted_gradient(G1);
  group2->set_weighted_gradient(G2);
}
