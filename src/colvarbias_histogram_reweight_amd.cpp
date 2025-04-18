// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include "colvarbias_histogram_reweight_amd.h"
#include "colvarproxy.h"
#include "colvars_memstream.h"

colvarbias_reweightaMD::colvarbias_reweightaMD(char const *key) : colvarbias_histogram(key) {}

colvarbias_reweightaMD::~colvarbias_reweightaMD() {}

int colvarbias_reweightaMD::init(std::string const &conf) {
  if (cvm::proxy->accelMD_enabled() == false) {
    cvm::error("Error: accelerated MD in your MD engine is not enabled.\n", COLVARS_INPUT_ERROR);
  }
  cvm::main()->cite_feature("reweightaMD colvar bias implementation (NAMD)");
  int baseclass_init_code = colvarbias_histogram::init(conf);
  get_keyval(conf, "CollectAfterSteps", start_after_steps, 0);
  get_keyval(conf, "CumulantExpansion", b_use_cumulant_expansion, true);
  get_keyval(conf, "WritePMFGradients", b_write_gradients, true);
  get_keyval(conf, "historyFreq", history_freq, 0);
  b_history_files = (history_freq > 0);
  grid_count.reset(new colvar_grid_scalar(colvars, nullptr, false, grid_conf));
  grid_count->request_actual_value();
  grid->request_actual_value();
  pmf_grid_exp_avg.reset(new colvar_grid_scalar(colvars, grid_count));
  if (b_write_gradients) {
    grad_grid_exp_avg.reset(new colvar_grid_gradient(colvars, nullptr, grid_count));
  }
  if (b_use_cumulant_expansion) {
    grid_dV.reset(new colvar_grid_scalar(colvars, grid_count));
    grid_dV_square.reset(new colvar_grid_scalar(colvars, grid_count));
    pmf_grid_cumulant.reset(new colvar_grid_scalar(colvars, grid_count));
    grid_dV->request_actual_value();
    grid_dV_square->request_actual_value();
    if (b_write_gradients) {
      grad_grid_cumulant.reset(new colvar_grid_gradient(colvars, nullptr, grid_count));
    }
  }
  previous_bin.assign(num_variables(), -1);
  return baseclass_init_code;
}

int colvarbias_reweightaMD::update() {
  colvarproxy *proxy = cvm::main()->proxy;
  int error_code = COLVARS_OK;
  if (cvm::step_relative() >= start_after_steps) {
    // update base class
    error_code |= colvarbias::update();

    if (cvm::debug()) {
      cvm::log("Updating histogram bias " + this->name);
    }

    if (cvm::step_relative() > 0) {
      previous_bin = bin;
    }

    // assign a valid bin size
    bin.assign(num_variables(), 0);

    if (colvar_array_size == 0) {
      // update indices for scalar values
      size_t i;
      for (i = 0; i < num_variables(); i++) {
        bin[i] = grid->current_bin_scalar(i);
      }

      if (grid->index_ok(previous_bin) && cvm::step_relative() > 0) {
        const cvm::real reweighting_factor = cvm::proxy->get_accelMD_factor();
        grid_count->acc_value(previous_bin, 1.0);
        grid->acc_value(previous_bin, reweighting_factor);
        if (b_use_cumulant_expansion) {
          const cvm::real dV = cvm::logn(reweighting_factor) *
            proxy->target_temperature() * proxy->boltzmann();
          grid_dV->acc_value(previous_bin, dV);
          grid_dV_square->acc_value(previous_bin, dV * dV);
        }
      }
    } else {
      // update indices for vector/array values
      size_t iv, i;
      for (iv = 0; iv < colvar_array_size; iv++) {
        for (i = 0; i < num_variables(); i++) {
          bin[i] = grid->current_bin_scalar(i, iv);
        }

      if (grid->index_ok(previous_bin) && cvm::step_relative() > 0) {
          const cvm::real reweighting_factor = cvm::proxy->get_accelMD_factor();
          grid_count->acc_value(previous_bin, 1.0);
          grid->acc_value(previous_bin, reweighting_factor);
          if (b_use_cumulant_expansion) {
            const cvm::real dV = cvm::logn(reweighting_factor) *
              proxy->target_temperature() * proxy->boltzmann();
            grid_dV->acc_value(previous_bin, dV);
            grid_dV_square->acc_value(previous_bin, dV * dV);
          }
        }
      }
    }
    previous_bin.assign(num_variables(), 0);

    error_code |= cvm::get_error();
  }
  return error_code;
}

int colvarbias_reweightaMD::write_output_files() {
  int error_code = COLVARS_OK;
  // error_code |= colvarbias_histogram::write_output_files();
  const std::string out_name_pmf = cvm::output_prefix() + "." +
                                   this->name + ".reweight";
  error_code |= write_exponential_reweighted_pmf(out_name_pmf);
  const std::string out_count_prefix = cvm::output_prefix() + "." +
                                       this->name;
  error_code |= write_count(out_count_prefix);
  const bool write_history = b_history_files &&
                             (cvm::step_absolute() % history_freq) == 0;
  if (write_history) {
    error_code |= write_exponential_reweighted_pmf(
      out_name_pmf + ".hist", true);
    error_code |= write_count(out_count_prefix + ".hist",
                              (cvm::step_relative() > 0));
  }
  if (b_use_cumulant_expansion) {
    const std::string out_name_cumulant_pmf = cvm::output_prefix() + "." +
                                              this->name + ".cumulant";
    error_code |= write_cumulant_expansion_pmf(out_name_cumulant_pmf);
    if (write_history) {
      error_code |= write_cumulant_expansion_pmf(
        out_name_cumulant_pmf + ".hist", true);
    }
  }
  error_code |= cvm::get_error();
  return error_code;
}

int colvarbias_reweightaMD::write_exponential_reweighted_pmf(
  const std::string& p_output_prefix, bool keep_open) {
  const std::string output_pmf = p_output_prefix + ".pmf";

  cvm::log("Writing the accelerated MD PMF file \"" + output_pmf + "\".\n");
  std::ostream &pmf_grid_os = cvm::proxy->output_stream(output_pmf, "PMF file");
  if (!pmf_grid_os) {
    return COLVARS_FILE_ERROR;
  }
  pmf_grid_exp_avg->copy_grid(*grid);
  // compute the average
  for (size_t i = 0; i < pmf_grid_exp_avg->raw_data_num(); ++i) {
    const double count = grid_count->value(i);
    if (count > 0) {
      const double tmp = pmf_grid_exp_avg->value(i);
      pmf_grid_exp_avg->set_value(i, tmp / count);
    }
  }
  hist_to_pmf(pmf_grid_exp_avg.get(), grid_count.get());
  pmf_grid_exp_avg->write_multicol(pmf_grid_os);
  if (!keep_open) {
    cvm::proxy->close_output_stream(output_pmf);
  }

  if (b_write_gradients) {
    const std::string output_grad = p_output_prefix + ".grad";
    cvm::log("Writing the accelerated MD gradients file \"" + output_grad +
             "\".\n");
    std::ostream &grad_grid_os = cvm::proxy->output_stream(output_grad, "gradient file");
    if (!grad_grid_os) {
      return COLVARS_FILE_ERROR;
    }
    for (std::vector<int> ix = grad_grid_exp_avg->new_index();
          grad_grid_exp_avg->index_ok(ix); grad_grid_exp_avg->incr(ix)) {
      for (size_t n = 0; n < grad_grid_exp_avg->multiplicity(); n++) {
        grad_grid_exp_avg->set_value(
          ix, pmf_grid_exp_avg->gradient_finite_diff(ix, n), n);
      }
    }
    grad_grid_exp_avg->write_multicol(grad_grid_os);
    if (!keep_open) {
      cvm::proxy->close_output_stream(output_grad);
    }
  }

  return COLVARS_OK;
}

int colvarbias_reweightaMD::write_cumulant_expansion_pmf(
  const std::string& p_output_prefix, bool keep_open) {
  const std::string output_pmf = p_output_prefix + ".pmf";
  cvm::log("Writing the accelerated MD PMF file using cumulant expansion: \"" + output_pmf + "\".\n");
  std::ostream &pmf_grid_cumulant_os = cvm::proxy->output_stream(output_pmf, "PMF file");
  if (!pmf_grid_cumulant_os) {
    return COLVARS_FILE_ERROR;
  }
  compute_cumulant_expansion_factor(grid_dV.get(), grid_dV_square.get(),
                                    grid_count.get(), pmf_grid_cumulant.get());
  hist_to_pmf(pmf_grid_cumulant.get(), grid_count.get());
  pmf_grid_cumulant->write_multicol(pmf_grid_cumulant_os);
  if (!keep_open) {
    cvm::proxy->close_output_stream(output_pmf);
  }

  if (b_write_gradients) {
    const std::string output_grad = p_output_prefix + ".grad";
    cvm::log("Writing the accelerated MD gradients file \"" + output_grad + "\".\n");
    std::ostream &grad_grid_os = cvm::proxy->output_stream(output_grad, "grad file");
    if (!grad_grid_os) {
      return COLVARS_FILE_ERROR;
    }
    for (std::vector<int> ix = grad_grid_cumulant->new_index();
          grad_grid_cumulant->index_ok(ix); grad_grid_cumulant->incr(ix)) {
      for (size_t n = 0; n < grad_grid_cumulant->multiplicity(); n++) {
        grad_grid_cumulant->set_value(
          ix, pmf_grid_cumulant->gradient_finite_diff(ix, n), n);
      }
    }
    grad_grid_cumulant->write_multicol(grad_grid_os);
    cvm::proxy->close_output_stream(output_grad);
  }
  return COLVARS_OK;
}

int colvarbias_reweightaMD::write_count(const std::string& p_output_prefix, bool keep_open) {
  const std::string output_name = p_output_prefix + ".count";
  cvm::log("Writing the accelerated MD count file \""+output_name+"\".\n");
  std::ostream &grid_count_os = cvm::proxy->output_stream(output_name, "count file");
  if (!grid_count_os) {
    return COLVARS_FILE_ERROR;
  }
  grid_count->write_multicol(grid_count_os);
  if (!keep_open) {
    cvm::proxy->close_output_stream(output_name);
  }
  return COLVARS_OK;
}

void colvarbias_reweightaMD::hist_to_pmf(
  colvar_grid_scalar* hist,
  const colvar_grid_scalar* hist_count) const
{
  colvarproxy *proxy = cvm::main()->proxy;
  if (hist->raw_data_num() == 0) return;
  const cvm::real kbt = proxy->boltzmann() * proxy->target_temperature();
  bool first_min_element = true;
  bool first_max_element = true;
  cvm::real min_element = 0;
  cvm::real max_element = 0;
  size_t i = 0;
  // the first loop: using logarithm to compute PMF
  for (i = 0; i < hist->raw_data_num(); ++i) {
    const cvm::real count = hist_count->value(i);
    if (count > 0) {
      const cvm::real x = hist->value(i);
      const cvm::real pmf_value = -1.0 * kbt * cvm::logn(x);
      hist->set_value(i, pmf_value);
      // find the minimum PMF value
      if (first_min_element) {
        // assign current PMF value to min_element at the first time
        min_element = pmf_value;
        first_min_element = false;
      } else {
        // if this is not the first time, then
        min_element = (pmf_value < min_element) ? pmf_value : min_element;
      }
      // do the same to the maximum
      if (first_max_element) {
        max_element = pmf_value;
        first_max_element = false;
      } else {
        max_element = (pmf_value > max_element) ? pmf_value : max_element;
      }
    }
  }
  // the second loop: bringing the minimum PMF value to zero
  for (i = 0; i < hist->raw_data_num(); ++i) {
    const cvm::real count = hist_count->value(i);
    if (count > 0) {
      // bins that have samples
      const cvm::real x = hist->value(i);
      hist->set_value(i, x - min_element);
    } else {
      hist->set_value(i, max_element - min_element);
    }
  }
}


void colvarbias_reweightaMD::compute_cumulant_expansion_factor(
  const colvar_grid_scalar* hist_dV,
  const colvar_grid_scalar* hist_dV_square,
  const colvar_grid_scalar* hist_count,
  colvar_grid_scalar* cumulant_expansion_factor) const
{
  colvarproxy *proxy = cvm::main()->proxy;
  const cvm::real beta = 1.0 / (proxy->boltzmann() * proxy->target_temperature());
  size_t i = 0;
  for (i = 0; i < hist_dV->raw_data_num(); ++i) {
    const cvm::real count = hist_count->value(i);
    if (count > 0) {
      const cvm::real dV_avg = hist_dV->value(i) / count;
      const cvm::real dV_square_avg = hist_dV_square->value(i) / count;
      const cvm::real factor = cvm::exp(beta * dV_avg + 0.5 * beta * beta * (dV_square_avg - dV_avg * dV_avg));
      cumulant_expansion_factor->set_value(i, factor);
    }
  }
}


template <typename OST> OST & colvarbias_reweightaMD::write_state_data_template_(OST& os)
{
  std::ios::fmtflags flags(os.flags());
  os.setf(std::ios::fmtflags(0), std::ios::floatfield);
  write_state_data_key(os, "grid");
  grid->write_raw(os, 8);
  write_state_data_key(os, "grid_count");
  grid_count->write_raw(os, 8);
  write_state_data_key(os, "grid_dV");
  grid_dV->write_raw(os, 8);
  write_state_data_key(os, "grid_dV_square");
  grid_dV_square->write_raw(os, 8);
  os.flags(flags);
  return os;
}


std::ostream & colvarbias_reweightaMD::write_state_data(std::ostream& os)
{
  return write_state_data_template_<std::ostream>(os);
}


cvm::memory_stream & colvarbias_reweightaMD::write_state_data(cvm::memory_stream& os)
{
  return write_state_data_template_<cvm::memory_stream>(os);
}


template <typename IST> IST & colvarbias_reweightaMD::read_state_data_template_(IST& is)
{
  if (! read_state_data_key(is, "grid")) {
    return is;
  }
  if (! grid->read_raw(is)) {
    return is;
  }
  if (! read_state_data_key(is, "grid_count")) {
    return is;
  }
  if (! grid_count->read_raw(is)) {
    return is;
  }
  if (! read_state_data_key(is, "grid_dV")) {
    return is;
  }
  if (! grid_dV->read_raw(is)) {
    return is;
  }
  if (! read_state_data_key(is, "grid_dV_square")) {
    return is;
  }
  if (! grid_dV_square->read_raw(is)) {
    return is;
  }
  return is;
}


std::istream & colvarbias_reweightaMD::read_state_data(std::istream& is)
{
  return read_state_data_template_<std::istream>(is);
}


cvm::memory_stream & colvarbias_reweightaMD::read_state_data(cvm::memory_stream& is)
{
  return read_state_data_template_<cvm::memory_stream>(is);
}
