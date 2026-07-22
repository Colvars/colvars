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

#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
#include "cuda/colvarcomp_coordnums_kernel.h"
#endif // defined (COLVARS_CUDA) || defined (COLVARS_HIP)

#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
class colvar::coordnum::coordnum_gpu_impl_t {
public:
  coordnum_gpu_impl_t(colvar::coordnum* cpu_coordnum_in):
    cvmodule(cpu_coordnum_in->cvmodule), cvc(cpu_coordnum_in),
    d_pairlist(nullptr), d_coordnum(nullptr),
    d_com_grad_tmp{nullptr, nullptr},
    d_com_grad_out{nullptr, nullptr},
    d_tbcount(nullptr), pairlist_transposed(false),
    numTiles(0), tileListsSize(0),
    d_tileLists(nullptr), d_tileListsStart(nullptr),
    d_tileListsLen(nullptr), h_coordnum(nullptr) {}
  ~coordnum_gpu_impl_t() {
    colvarproxy* p = cvmodule->proxy;
    p->deallocate_device(&d_pairlist);
    p->deallocate_device(&d_coordnum);
    p->deallocate_device(&d_tbcount);
    p->deallocate_host(&h_coordnum);
    p->deallocate_device(&d_com_grad_tmp[0]);
    p->deallocate_device(&d_com_grad_tmp[1]);
    p->deallocate_device(&d_com_grad_out[0]);
    p->deallocate_device(&d_com_grad_out[1]);
    p->deallocate_device(&d_tileLists);
    p->deallocate_device(&d_tileListsStart);
    p->deallocate_device(&d_tileListsLen);
  }
  int init() {
    int error_code = COLVARS_OK;
    colvarproxy* p = cvmodule->proxy;
    error_code |= p->reallocate_device(&d_coordnum, 1);
    error_code |= p->reallocate_device(&d_tbcount, 1);
    error_code |= p->reallocate_host(&h_coordnum, 1);
    if (cvc->pairlist != nullptr && cvc->num_pairs > 0) {
      error_code |= p->reallocate_device(&d_pairlist, cvc->num_pairs);
      error_code |= p->clear_device_array(d_pairlist, cvc->num_pairs);
    }
    error_code |= p->clear_device_array(d_coordnum, 1);
    error_code |= p->clear_device_array(d_tbcount, 1);
    h_coordnum[0] = 0;
    if (cvc->function_type() != "selfCoordNum") {
      if (cvc->b_group1_center_only || cvc->b_group2_center_only) {
        error_code |= p->reallocate_device(&d_com_grad_tmp[0], 1);
        error_code |= p->reallocate_device(&d_com_grad_tmp[1], 1);
        error_code |= p->reallocate_device(&d_com_grad_out[0], 1);
        error_code |= p->reallocate_device(&d_com_grad_out[1], 1);
        error_code |= p->clear_device_array(d_com_grad_tmp[0], 1);
        error_code |= p->clear_device_array(d_com_grad_tmp[1], 1);
        error_code |= p->clear_device_array(d_com_grad_out[0], 1);
        error_code |= p->clear_device_array(d_com_grad_out[1], 1);
      }
      pairlist_transposed = cvc->group1->size() > cvc->group2->size();
    }
    return error_code;
  }
  int calc_value_two_groups(int flags) {
    int error_code = COLVARS_OK;
    error_code |= checkGPUError(cudaStreamWaitEvent(
      cvc->get_stream(), cvmodule->proxy->get_event(colvarproxy_gpu::event_type::update_lattice)));
    error_code |= checkGPUError(cudaStreamWaitEvent(
      cvc->get_stream(), cvc->group1->get_gpu_atom_group()->get_event(
        colvars_gpu::colvaratoms_gpu::event_type::read_and_calculate)));
    error_code |= checkGPUError(cudaStreamWaitEvent(
      cvc->get_stream(), cvc->group2->get_gpu_atom_group()->get_event(
        colvars_gpu::colvaratoms_gpu::event_type::read_and_calculate)));
    error_code |= colvars_gpu::calc_value_coordnum_two_groups(
      cvc->group1->get_gpu_atom_group()->get_gpu_buffers().d_atoms_pos,
      cvc->group2->get_gpu_atom_group()->get_gpu_buffers().d_atoms_pos,
      cvc->group1->size(), cvc->group2->size(), cvc->en, cvc->ed,
      cvc->inv_r0_vec, cvc->inv_r0sq_vec,
#if 0
      // TODO: Wait for https://github.com/Colvars/colvars/pull/919
#endif
      cvc->group1->get_gpu_atom_group()->get_gpu_buffers().d_atoms_grad,
      cvc->group2->get_gpu_atom_group()->get_gpu_buffers().d_atoms_grad,
      cvc->tolerance, cvc->tolerance_l2_max, d_pairlist, d_tbcount,
      d_coordnum, h_coordnum, flags, cvc->get_stream(), cvmodule);
    error_code |= checkGPUError(cudaEventRecord(
      cvc->get_event(cvc::event_type::calc_value), cvc->get_stream()));
    if (flags & colvar::coordnum::ef_gradients) {
      error_code |= checkGPUError(cudaEventRecord(
      cvc->get_event(cvc::event_type::calc_gradients), cvc->get_stream()));
    }
    return error_code;
  }
  int calc_value_group_to_com(int flags) {
    // Requires b_group1_center_only != b_group2_center_only
    if (cvc->b_group1_center_only == cvc->b_group2_center_only) {
      return cvmodule->error("calc_value_group_to_com is called incorrectly.\n",
                             COLVARS_BUG_ERROR);
    }
    int error_code = COLVARS_OK;
    auto group = cvc->b_group1_center_only ? cvc->group2 : cvc->group1;
    auto group_com = cvc->b_group1_center_only ? cvc->group1 : cvc->group2;
    error_code |= checkGPUError(cudaStreamWaitEvent(
      cvc->get_stream(), cvmodule->proxy->get_event(colvarproxy_gpu::event_type::update_lattice)));
    error_code |= checkGPUError(cudaStreamWaitEvent(
      cvc->get_stream(), cvc->group1->get_gpu_atom_group()->get_event(
        colvars_gpu::colvaratoms_gpu::event_type::read_and_calculate)));
    error_code |= checkGPUError(cudaStreamWaitEvent(
      cvc->get_stream(), cvc->group2->get_gpu_atom_group()->get_event(
        colvars_gpu::colvaratoms_gpu::event_type::read_and_calculate)));
    error_code |= colvars_gpu::calc_value_coordnum_group_to_com(
      group->get_gpu_atom_group()->get_gpu_buffers().d_atoms_pos,
      group_com->get_gpu_atom_group()->get_gpu_buffers().d_com,
      cvc->group1->size(), cvc->en, cvc->ed, cvc->inv_r0_vec,
      cvc->inv_r0sq_vec,
      group->get_gpu_atom_group()->get_gpu_buffers().d_atoms_grad,
      cvc->tolerance, cvc->tolerance_l2_max, d_pairlist, d_tbcount,
      d_com_grad_tmp[0], d_com_grad_out[0], d_coordnum, h_coordnum,
      flags, cvc->get_stream(), cvmodule);
    error_code |= checkGPUError(cudaEventRecord(
      cvc->get_event(cvc::event_type::calc_value), cvc->get_stream()));
    if (flags & colvar::coordnum::ef_gradients) {
      // Set weighted gradients
      error_code |= colvars_gpu::colvaratoms_gpu::set_weighted_gradient_gpu(
        group_com, d_com_grad_out[0], cvc->get_stream());
      error_code |= checkGPUError(cudaEventRecord(
        cvc->get_event(cvc::event_type::calc_gradients), cvc->get_stream()));
    }
    return error_code;
  }
  int calc_value_two_coms(int flags) {
    const bool check = (cvc->b_group1_center_only == true) && (cvc->b_group2_center_only == true);
    if (!check) {
      return cvmodule->error("calc_value_two_coms is called incorrectly.\n",
                             COLVARS_BUG_ERROR);
    }
    int error_code = COLVARS_OK;
    error_code |= checkGPUError(cudaStreamWaitEvent(
      cvc->get_stream(), cvmodule->proxy->get_event(colvarproxy_gpu::event_type::update_lattice)));
    error_code |= checkGPUError(cudaStreamWaitEvent(
      cvc->get_stream(), cvc->group1->get_gpu_atom_group()->get_event(
        colvars_gpu::colvaratoms_gpu::event_type::read_and_calculate)));
    error_code |= checkGPUError(cudaStreamWaitEvent(
      cvc->get_stream(), cvc->group2->get_gpu_atom_group()->get_event(
        colvars_gpu::colvaratoms_gpu::event_type::read_and_calculate)));
    error_code |= colvars_gpu::calc_value_coordnum_com_to_com(
      cvc->group1->get_gpu_atom_group()->get_gpu_buffers().d_com,
      cvc->group2->get_gpu_atom_group()->get_gpu_buffers().d_com,
      cvc->en, cvc->ed, cvc->inv_r0_vec, cvc->inv_r0sq_vec,
      cvc->tolerance, cvc->tolerance_l2_max, d_pairlist,
      d_com_grad_out[0], d_com_grad_out[1], h_coordnum, flags, cvc->get_stream(), cvmodule);
    error_code |= checkGPUError(cudaEventRecord(
      cvc->get_event(cvc::event_type::calc_value), cvc->get_stream()));
    if (flags & colvar::coordnum::ef_gradients) {
      // Set weighted gradients
      error_code |= colvars_gpu::colvaratoms_gpu::set_weighted_gradient_gpu(
        cvc->group1, d_com_grad_out[0], cvc->get_stream());
      error_code |= colvars_gpu::colvaratoms_gpu::set_weighted_gradient_gpu(
        cvc->group2, d_com_grad_out[1], cvc->get_stream());
      error_code |= checkGPUError(cudaEventRecord(
        cvc->get_event(cvc::event_type::calc_gradients), cvc->get_stream()));
    }
    return error_code;
  }
  // NOTE: Maybe it is better to build tile lists on GPU if we
  // want to exclude non-interacting tiles periodically like NAMD.
  int buildTileLists(const unsigned int tileSize, cudaStream_t stream) {
    int error_code = COLVARS_OK;
    colvarproxy* p = cvmodule->proxy;
    const unsigned int numAtoms = cvc->group1->size();
    numTiles = (numAtoms + tileSize - 1) / tileSize;
    const unsigned int numTileInteractions = numTiles * (numTiles - 1) / 2;
    const unsigned int maxNumInteractionsPerTile =
      (numTileInteractions + numTiles - 1) / numTiles;
    tileListsSize = numTiles * maxNumInteractionsPerTile;
    tileListsLen.resize(numTiles, 0);
    tileListsStart.resize(numTiles, 0);
    tileLists.resize(tileListsSize);
    /**
      for (unsigned int i = 0; i < numTiles; ++i) {
        tileListsStart[i] = i * maxNumInteractionsPerTile;
        const unsigned int jStart = i + 1;
        const unsigned int jEnd = std::min(numTiles, jStart + maxNumInteractionsPerTile);
        for (unsigned int j = jStart; j < jEnd; ++j) {
          const unsigned int offset = tileListsLen[i]++;
          const unsigned int pos = i * maxNumInteractionsPerTile;
          tileLists[pos+offset] = j;
        }
        for (unsigned int j = jEnd; j < numTiles; ++j) {
          const unsigned int offset = tileListsLen[j]++;
          const unsigned int pos = j * maxNumInteractionsPerTile;
          tileLists[pos+offset] = i;
        }
      }
      * @note The following loop is a more refined version of the loop above
    */
    for (unsigned int i = 0; i < numTiles; ++i) {
      const unsigned int jStart = i + 1;
      const unsigned int jEnd = std::min(numTiles, jStart + maxNumInteractionsPerTile);
      const unsigned int posStart = i * maxNumInteractionsPerTile;
      tileListsLen[i] += jEnd - jStart;
      tileListsStart[i] = posStart;
      for (unsigned int j = jStart; j < jEnd; ++j) {
        const unsigned int offset = j - jStart;
        tileLists[posStart+offset] = j;
      }
      if (i > maxNumInteractionsPerTile) {
        tileListsLen[i] += i - maxNumInteractionsPerTile;
        for (unsigned int k = 0; k < i - maxNumInteractionsPerTile; ++k) {
          const unsigned int offset = jEnd - jStart + k;
          tileLists[posStart+offset] = k;
        }
      }
    }
#if 0
    if constexpr (cvm::debug()) {
      for (unsigned int i = 0; i < numTiles; ++i) {
        const unsigned int pos = tileListsStart[i];
        const unsigned int len = tileListsLen[i];
        const auto begin = tileLists.begin() + pos;
        const auto end = begin + len;
        std::string s = "Tilelist[" + cvm::to_str(int(i)) + "], start = " +
                        cvm::to_str(int(pos)) + " size = " + cvm::to_str(int(len)) +
                        ", data = [";
        for (auto it = begin; it != end; ++it) {
          const unsigned int tile = *it;
          s += " " + cvm::to_str((int)tile);
        }
        s += "]\n";
        cvmodule->log(s);
      }
    }
#endif
    error_code |= p->deallocate_device_async(&d_tileLists, stream);
    error_code |= p->deallocate_device_async(&d_tileListsStart, stream);
    error_code |= p->deallocate_device_async(&d_tileListsLen, stream);
    error_code |= p->allocate_device_async(&d_tileLists, tileListsSize, stream);
    error_code |= p->allocate_device_async(&d_tileListsStart, numTiles, stream);
    error_code |= p->allocate_device_async(&d_tileListsLen, numTiles, stream);
    error_code |= p->copy_HtoD_async(tileLists.data(), d_tileLists, tileListsSize, stream);
    error_code |= p->copy_HtoD_async(tileListsStart.data(), d_tileListsStart, numTiles, stream);
    error_code |= p->copy_HtoD_async(tileListsLen.data(), d_tileListsLen, numTiles, stream);
    return error_code;
  }
  int calc_value_self_group(int flags) {
    int error_code = COLVARS_OK;
    if (d_tileLists == nullptr) {
      colvarproxy* p = cvmodule->proxy;
      const unsigned int gpu_warp_size = p->gpu_warp_size();
      error_code |= buildTileLists(gpu_warp_size, cvc->get_stream());
    }
    error_code |= checkGPUError(cudaStreamWaitEvent(
      cvc->get_stream(), cvmodule->proxy->get_event(colvarproxy_gpu::event_type::update_lattice)));
    error_code |= checkGPUError(cudaStreamWaitEvent(
      cvc->get_stream(), cvc->group1->get_gpu_atom_group()->get_event(
        colvars_gpu::colvaratoms_gpu::event_type::read_and_calculate)));
    error_code |= colvars_gpu::calc_value_coordnum_self_group(
      cvc->group1->get_gpu_atom_group()->get_gpu_buffers().d_atoms_pos,
      cvc->group1->size(), cvc->en, cvc->ed, cvc->inv_r0_vec, cvc->inv_r0sq_vec,
      cvc->group1->get_gpu_atom_group()->get_gpu_buffers().d_atoms_grad,
      d_tileLists, d_tileListsStart, d_tileListsLen,
      cvc->tolerance, cvc->tolerance_l2_max, d_pairlist,
      d_tbcount, d_coordnum, h_coordnum, flags, cvc->get_stream(), cvmodule);
    error_code |= checkGPUError(cudaEventRecord(
      cvc->get_event(cvc::event_type::calc_value), cvc->get_stream()));
    if (flags & colvar::coordnum::ef_gradients) {
      error_code |= checkGPUError(cudaEventRecord(
        cvc->get_event(cvc::event_type::calc_gradients), cvc->get_stream()));
    }
    return error_code;
  }
private:
  colvarmodule* cvmodule;
  colvar::coordnum* cvc;
  bool* d_pairlist;
  cvm::real* d_coordnum;
  cvm::rvector* d_com_grad_tmp[2];
  cvm::rvector* d_com_grad_out[2];
  unsigned int* d_tbcount;
  bool pairlist_transposed;
  /**
   * @name Tile lists for self coordnum
   */
  ///@[
  using tl_vec_t = std::vector<unsigned int, colvars_gpu::CudaHostAllocator<unsigned int>>;
  tl_vec_t tileListsLen;
  tl_vec_t tileListsStart;
  tl_vec_t tileLists;
  unsigned int numTiles;
  unsigned int tileListsSize;
  unsigned int* d_tileLists;
  unsigned int* d_tileListsStart;
  unsigned int* d_tileListsLen;
  ///@]
public:
  cvm::real* h_coordnum;
};
#endif


colvar::coordnum::coordnum()
{
  set_function_type("coordNum");
  x.type(colvarvalue::type_scalar);
  cvm::real const r0 = cvmodule->proxy->angstrom_to_internal(4.0);
  update_cutoffs({r0, r0, r0});
  b_use_internal_pbc = cvm::main()->proxy->use_internal_pbc();
  // Boundaries will be set later, when the number of pairs is known
}


void colvar::coordnum::update_cutoffs(cvm::rvector const &r0_vec_i)
{
  r0_vec = r0_vec_i;

  inv_r0_vec = {
    1.0 / r0_vec.x,
    1.0 / r0_vec.y,
    1.0 / r0_vec.z
  };

  inv_r0sq_vec = {
    inv_r0_vec.x * inv_r0_vec.x,
    inv_r0_vec.y * inv_r0_vec.y,
    inv_r0_vec.z * inv_r0_vec.z
  };
}


int colvar::coordnum::init(std::string const &conf)
{
  int error_code = cvc::init(conf);

  get_keyval(conf, "useInternalPBC", b_use_internal_pbc, b_use_internal_pbc);

  group1 = parse_group(conf, "group1");

  if (!group1) {
    return error_code | COLVARS_INPUT_ERROR;
  }

  if (group1->b_dummy) {
    error_code |= cvmodule->error("Error: group1 may not be a dummy atom\n", COLVARS_INPUT_ERROR);
  }

  if (function_type() != "selfCoordNum") {

    group2 = parse_group(conf, "group2");
    if (!group2) {
      return error_code | COLVARS_INPUT_ERROR;
    }

    if (int atom_number = cvm::atom_group::overlap(*group1, *group2)) {
      error_code |= cvmodule->error("Error: group1 and group2 share a common atom (number: " +
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
          cvmodule->error("Error: cannot specify \"cutoff\" and \"cutoff3\" at the same time.\n",
                          COLVARS_INPUT_ERROR);
    }

    // remove meaningless negative signs
    if (r0_vec.x < 0.0) r0_vec.x *= -1.0;
    if (r0_vec.y < 0.0) r0_vec.y *= -1.0;
    if (r0_vec.z < 0.0) r0_vec.z *= -1.0;

    update_cutoffs(r0_vec);

  } else {
    if (b_redefined_cutoff) {
      update_cutoffs({r0, r0, r0});
    }
  }

  get_keyval(conf, "expNumer", en, en);
  get_keyval(conf, "expDenom", ed, ed);

  if ( (en%2) || (ed%2) ) {
    error_code |= cvmodule->error("Error: odd exponent(s) provided, can only use even ones.\n",
                             COLVARS_INPUT_ERROR);
  }

  if ( (en <= 0) || (ed <= 0) ) {
    error_code |= cvmodule->error("Error: negative exponent(s) provided.\n",
                             COLVARS_INPUT_ERROR);
  }

  if (!is_enabled(f_cvc_pbc_minimum_image)) {
    cvmodule->log("Warning: only minimum-image distances are used by this variable.\n");
  }

  if (function_type() != "groupCoord") {
    // All coordNum variables may benefit from a pairlist, except groupCoord
    get_keyval(conf, "tolerance", tolerance, tolerance);
    if (tolerance > 0) {
      cvmodule->cite_feature("coordNum pairlist");
      compute_tolerance_l2_max();
      get_keyval(conf, "pairListFrequency", pairlist_freq, pairlist_freq);
      if ( ! (pairlist_freq > 0) ) {
        return cvmodule->error("Error: non-positive pairlistfrequency provided.\n",
                               COLVARS_INPUT_ERROR);
        // return and do not allocate the pairlists below
      }
      pairlist.reset(new bool[num_pairs]);
      auto *pairlist_elem = pairlist.get();
      for (size_t ip = 0; ip < num_pairs; ip++, pairlist_elem++) {
        *pairlist_elem = true;
      }
    }
  }

  if (cvmodule->proxy->get_smp_mode() == colvarproxy_smp::smp_mode_t::gpu) {
#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
    coordnum_gpu_impl = std::unique_ptr<coordnum_gpu_impl_t>(new coordnum_gpu_impl_t(this));
    error_code |= coordnum_gpu_impl->init();
#endif
  }
  return error_code;
}


colvar::coordnum::~coordnum() {}


void colvar::coordnum::compute_tolerance_l2_max()
{
  cvm::real l2 = 1.001;
  cvm::real F = 0.0;
  cvm::real dFdl2 = 0.0;
  constexpr size_t num_iters_max = 1000000;
  constexpr cvm::real result_tol = 1.0e-6;
  constexpr cvm::real dF_tol = 1.0e-9;
  size_t i;
  // Find the value of l2 such that F(l2) = 0 using the Newton method
  for (i = 0; i < num_iters_max; i++) {
    F = switching_function<ef_use_pairlist | ef_gradients>(l2, dFdl2, en, ed, tolerance);
    if ((std::fabs(F) < result_tol) || (std::fabs(dFdl2) < dF_tol)) {
      break;
    }
    l2 -= F / dFdl2;
  }
  tolerance_l2_max = l2;
  if (cvm::debug()) {
    cvmodule->log("Found max valid l2 in " + cvm::to_str(i + 1) +
                  " iterations, result = " + cvm::to_str(l2) + " f(result) = " + cvm::to_str(F));
  }
}


template <bool use_group1_com, bool use_group2_com, int flags>
void inline colvar::coordnum::main_loop()
{
  size_t const group1_num_coords = use_group1_com ? 1 : group1->size();
  size_t const group2_num_coords = use_group2_com ? 1 : group2->size();

  cvm::atom_pos const group1_com = group1->center_of_mass();
  cvm::atom_pos const group2_com = group2->center_of_mass();
  cvm::rvector group1_com_grad, group2_com_grad;

  bool *pairlist_elem = pairlist.get();

  for (size_t i = 0; i < group1_num_coords; ++i) {

    cvm::real const x1 = use_group1_com ? group1_com.x : group1->pos_x(i);
    cvm::real const y1 = use_group1_com ? group1_com.y : group1->pos_y(i);
    cvm::real const z1 = use_group1_com ? group1_com.z : group1->pos_z(i);

    cvm::real &gx1 = use_group1_com ? group1_com_grad.x : group1->grad_x(i);
    cvm::real &gy1 = use_group1_com ? group1_com_grad.y : group1->grad_y(i);
    cvm::real &gz1 = use_group1_com ? group1_com_grad.z : group1->grad_z(i);

    for (size_t j = 0; j < group2_num_coords; ++j) {

      cvm::real const x2 = use_group2_com ? group2_com.x : group2->pos_x(j);
      cvm::real const y2 = use_group2_com ? group2_com.y : group2->pos_y(j);
      cvm::real const z2 = use_group2_com ? group2_com.z : group2->pos_z(j);

      cvm::real &gx2 = use_group2_com ? group2_com_grad.x : group2->grad_x(j);
      cvm::real &gy2 = use_group2_com ? group2_com_grad.y : group2->grad_y(j);
      cvm::real &gz2 = use_group2_com ? group2_com_grad.z : group2->grad_z(j);

      bool const within =
          ((flags & ef_use_pairlist) && (*pairlist_elem || (flags & ef_rebuild_pairlist))) ||
          !(flags & ef_use_pairlist);

      cvm::real const partial = within ?
        (b_use_internal_pbc ?
         compute_pair_coordnum<flags | ef_use_internal_pbc>(inv_r0_vec, inv_r0sq_vec, en, ed,
                                                            x1, y1, z1, x2, y2, z2,
                                                            gx1, gy1, gz1, gx2, gy2, gz2,
                                                            tolerance, tolerance_l2_max,
                                                            cvmodule) :
         compute_pair_coordnum<flags>(inv_r0_vec, inv_r0sq_vec, en, ed,
                                      x1, y1, z1, x2, y2, z2,
                                      gx1, gy1, gz1, gx2, gy2, gz2,
                                      tolerance, tolerance_l2_max,
                                      cvmodule) ) :
        0.0;

      if ((flags & ef_use_pairlist) && (flags & ef_rebuild_pairlist)) {
        *pairlist_elem = partial > 0.0 ? true : false;
      }

      x.real_value += partial;

      if (flags & ef_use_pairlist) {
        pairlist_elem++;
      }
    }
  }

  if (use_group1_com) {
    group1->set_weighted_gradient(group1_com_grad);
  }
  if (use_group2_com) {
    group2->set_weighted_gradient(group2_com_grad);
  }
}


template <bool use_group1_com, bool use_group2_com, int compute_flags>
int colvar::coordnum::compute_coordnum()
{
  bool const use_pairlist = pairlist.get();
  bool const rebuild_pairlist = use_pairlist && (cvmodule->step_relative() % pairlist_freq == 0);

  if (use_pairlist) {
    if (rebuild_pairlist) {
      constexpr int flags = compute_flags | ef_use_pairlist | ef_rebuild_pairlist;
      main_loop<use_group1_com, use_group2_com, flags>();
    } else {
      constexpr int flags = compute_flags | ef_use_pairlist;
      main_loop<use_group1_com, use_group2_com, flags>();
    }
  } else {
    constexpr int flags = compute_flags;
    main_loop<use_group1_com, use_group2_com, flags>();
  }

  return COLVARS_OK;
}


void colvar::coordnum::calc_value()
{
  x.real_value = 0.0;
  if (is_enabled(f_cvc_gradient)) {

    constexpr int flags = ef_gradients;

    if (b_group1_center_only) {
      if (b_group2_center_only) {
        compute_coordnum<true, true, flags>();
      } else {
        compute_coordnum<true, false, flags>();
      }
    } else {
      if (b_group2_center_only) {
        compute_coordnum<false, true, flags>();
      } else {
        compute_coordnum<false, false, flags>();
      }
    }

  } else {

    constexpr int flags = ef_null;

    if (b_group1_center_only) {
      if (b_group2_center_only) {
        compute_coordnum<true, true, flags>();
      } else {
        compute_coordnum<true, false, flags>();
      }
    } else {
      if (b_group2_center_only) {
        compute_coordnum<false, true, flags>();
      } else {
        compute_coordnum<false, false, flags>();
      }
    }
  }
}


void colvar::coordnum::calc_gradients()
{
  // Gradients are computed by calc_value() if f_cvc_gradients is enabled
}

#if defined (COLVARS_CUDA) || defined (COLVARS_HIP)
int colvar::coordnum::calc_value_gpu() {
  int error_code = COLVARS_OK;
  const bool use_pairlist = pairlist.get();
  const bool rebuild_pairlist = use_pairlist && (cvmodule->step_relative() % pairlist_freq == 0);
  const bool gradients = is_enabled(f_cvc_gradient);
  const int flags = (use_pairlist ? ef_use_pairlist : 0) +
                    (rebuild_pairlist ? ef_rebuild_pairlist : 0) +
                    (gradients ? ef_gradients : 0);
  if (function_type() != "selfCoordNum") {
    if (!b_group1_center_only && !b_group2_center_only) {
      // Compute coordnum between two different sets of atoms
      error_code |= coordnum_gpu_impl->calc_value_two_groups(flags);
    } else {
      if (b_group1_center_only != b_group2_center_only) {
        // Compute coordnum between a group of atoms and the COM of another group
        error_code |= coordnum_gpu_impl->calc_value_group_to_com(flags);
      } else {
        // Compute coordnum between two COMs
        error_code |= coordnum_gpu_impl->calc_value_two_coms(flags);
      }
    }
  } else {
    // Compute coordnum within a group
    error_code |= coordnum_gpu_impl->calc_value_self_group(flags);
  }
  return error_code;
}

int colvar::coordnum::calc_value_after_gpu() {
  int error_code = COLVARS_OK;
  error_code |= checkGPUError(cudaEventSynchronize(get_event(cvc::event_type::calc_value)));
  x.real_value = *coordnum_gpu_impl->h_coordnum;
  return error_code;
}
#endif

// h_bond member functions

colvar::h_bond::h_bond()
{
  cvm::real const r0 = cvmodule->proxy->angstrom_to_internal(3.3);
  r0_vec = {r0, r0, r0};
  set_function_type("hBond");
  x.type(colvarvalue::type_scalar);
  init_scalar_boundaries(0.0, 1.0);
}


int colvar::h_bond::init(std::string const &conf)
{
  int error_code = cvc::init(conf);

  if (cvm::debug())
    cvmodule->log("Initializing h_bond object.\n");

  set_function_type("hBond");
  x.type(colvarvalue::type_scalar);
  init_scalar_boundaries(0.0, 1.0);

  int a_num = -1, d_num = -1;
  get_keyval(conf, "acceptor", a_num, a_num);
  get_keyval(conf, "donor",    d_num, a_num);

  if ( (a_num == -1) || (d_num == -1) ) {
    error_code |= cvmodule->error("Error: either acceptor or donor undefined.\n", COLVARS_INPUT_ERROR);
  }

  register_atom_group(new cvm::atom_group);
  {
    colvarproxy* const p = cvmodule->proxy;
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
    error_code |= cvmodule->error("Error: odd exponent(s) provided, can only use even ones.\n",
                             COLVARS_INPUT_ERROR);
  }

  if ((en <= 0) || (ed <= 0)) {
    error_code |= cvmodule->error("Error: negative exponent(s) provided.\n", COLVARS_INPUT_ERROR);
  }

  if (cvm::debug())
    cvmodule->log("Done initializing h_bond object.\n");

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
  constexpr int flags = coordnum::ef_null;
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

  x.real_value = coordnum::compute_pair_coordnum<flags>(inv_r0_vec, inv_r0sq_vec, en, ed,
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
                                                        0.0, 1.0e20,
                                                        cvmodule);
  // Skip the gradient
}


void colvar::h_bond::calc_gradients()
{
  int constexpr flags = coordnum::ef_gradients;
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
  coordnum::compute_pair_coordnum<flags>(inv_r0_vec, inv_r0sq_vec, en, ed,
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
                                         0.0, 1.0e20,
                                         cvmodule);
}


colvar::selfcoordnum::selfcoordnum()
{
  set_function_type("selfCoordNum");
}


template <int flags> inline void colvar::selfcoordnum::selfcoordnum_sequential_loop()
{
  size_t const n = group1->size();
  bool *pairlist_elem = pairlist.get();

  for (size_t i = 0; i < n - 1; i++) {

    cvm::real const x1 = group1->pos_x(i);
    cvm::real const y1 = group1->pos_y(i);
    cvm::real const z1 = group1->pos_z(i);
    cvm::real &gx1 = group1->grad_x(i);
    cvm::real &gy1 = group1->grad_y(i);
    cvm::real &gz1 = group1->grad_z(i);

    for (size_t j = i + 1; j < n; j++) {

      cvm::real const x2 = group1->pos_x(j);
      cvm::real const y2 = group1->pos_y(j);
      cvm::real const z2 = group1->pos_z(j);
      cvm::real &gx2 = group1->grad_x(j);
      cvm::real &gy2 = group1->grad_y(j);
      cvm::real &gz2 = group1->grad_z(j);

      bool const within =
        ((flags & ef_use_pairlist) && (*pairlist_elem || (flags & ef_rebuild_pairlist))) ||
        !(flags & ef_use_pairlist);

      cvm::real const partial = within ?
        (b_use_internal_pbc ?
         compute_pair_coordnum<flags | ef_use_internal_pbc>(inv_r0_vec, inv_r0sq_vec, en, ed,
                                                            x1, y1, z1, x2, y2, z2,
                                                            gx1, gy1, gz1, gx2, gy2, gz2,
                                                            tolerance, tolerance_l2_max,
                                                            cvmodule) :
         compute_pair_coordnum<flags>(inv_r0_vec, inv_r0sq_vec, en, ed,
                                      x1, y1, z1, x2, y2, z2,
                                      gx1, gy1, gz1, gx2, gy2, gz2,
                                      tolerance, tolerance_l2_max,
                                      cvmodule) ) :
        0.0;

      if ((flags & ef_use_pairlist) && (flags & ef_rebuild_pairlist)) {
        *pairlist_elem = partial > 0.0 ? true : false;
      }

      x.real_value += partial;

      if (flags & ef_use_pairlist) {
        pairlist_elem++;
      }
    }
  }
}


template<int compute_flags> int colvar::selfcoordnum::compute_selfcoordnum()
{
  bool const use_pairlist = pairlist.get();
  bool const rebuild_pairlist = use_pairlist && (cvmodule->step_relative() % pairlist_freq == 0);

  if (use_pairlist) {
    if (rebuild_pairlist) {
      int constexpr flags = compute_flags | ef_use_pairlist | ef_rebuild_pairlist;
      selfcoordnum_sequential_loop<flags>();
    } else {
      int constexpr flags = compute_flags | ef_use_pairlist;
      selfcoordnum_sequential_loop<flags>();
    }
  } else {
    int constexpr flags = compute_flags | ef_null;
    selfcoordnum_sequential_loop<flags>();
  }
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


colvar::groupcoordnum::groupcoordnum() { set_function_type("groupCoord"); }


void colvar::groupcoordnum::calc_value()
{
  x.real_value = 0.0;
  if (is_enabled(f_cvc_gradient)) {
    constexpr int flags = ef_gradients;
    compute_coordnum<true, true, flags>();
  } else {
    constexpr int flags = ef_null;
    compute_coordnum<true, true, flags>();
  }
}


void colvar::groupcoordnum::calc_gradients()
{
  // Gradients are computed by calc_value() if f_cvc_gradients is enabled
}
