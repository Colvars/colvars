#ifdef TORCH

#include "colvarmodule.h"
#include "colvar.h"
#include "colvarvalue.h"
#include "colvarparse.h"
#include "colvarcomp.h"

colvar::torchANN::torchANN(std::string const &conf): linearCombination(conf) {
  set_function_type("torchANN");

  x.type(colvarvalue::type_scalar);
  enable(f_cvc_scalar);

  if (period != 0.0) {
    enable(f_cvc_periodic);
  }

  if ((wrap_center != 0.0) && !is_enabled(f_cvc_periodic)) {
    cvm::error("Error: wrapAround was defined in a torchANN component,"
                " but its period has not been set.\n");
    return;
  }

  std::string model_file ;
  get_keyval(conf, "model_file", model_file, std::string(""));
  get_keyval(conf, "m_output_index", m_output_index, 0);

  try {
    nn = torch::jit::load(model_file);
    cvm::log("model loaded.") ;
  } catch (...) {
    cvm::error("Error: couldn't load libtorch model.\n");
  }

  cvc_indices.resize(cv.size(),0);

  size_t num_inputs = 0;
  // compute total number of inputs of neural network 
  for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) 
  {
      num_inputs += cv[i_cv]->value().size() ;
      if (i_cv < cv.size() - 1) 
	cvc_indices[i_cv+1] = num_inputs;
  }

  // initialize the input tensor 
  input_tensor = torch::zeros({1,(long int) num_inputs}, torch::TensorOptions().dtype(torch::kFloat32).requires_grad(true));
}

colvar::torchANN::~torchANN() {
}

void colvar::torchANN::calc_value() {

  for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) 
      cv[i_cv]->calc_value();

  // set input tensor with no_grad 
  {
    torch::NoGradGuard no_grad;
    size_t l = 0;
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
	const colvarvalue& current_cv_value = cv[i_cv]->value();
	if (current_cv_value.type() == colvarvalue::type_scalar) {
	    input_tensor[0][l++] = cv[i_cv]->sup_coeff * (cvm::pow(current_cv_value.real_value, cv[i_cv]->sup_np));
	} else {  
	    for (size_t j_elem = 0; j_elem < current_cv_value.size(); ++j_elem) 
		input_tensor[0][l++] = cv[i_cv]->sup_coeff * current_cv_value[j_elem];
	}
    }
  }
  if (input_tensor.grad().defined())
    input_tensor.grad().zero_();

  std::vector<torch::jit::IValue> inputs={input_tensor};

  // evaluate the value of function
  nn_outputs = nn.forward(inputs).toTensor()[0][m_output_index];

  nn_outputs.backward({}, false, false);
  input_grad = input_tensor.grad()[0];

  x = nn_outputs.item<double>() ;

  this->wrap(x);
}

void colvar::torchANN::calc_gradients() {

  for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
      cv[i_cv]->calc_gradients();
      if (cv[i_cv]->is_enabled(f_cvc_explicit_gradient)) {
	  const cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
	  // get the initial index of this cvc
	  size_t l = cvc_indices[i_cv];
	  for (size_t j_elem = 0; j_elem < cv[i_cv]->value().size(); ++j_elem) {
	      // get derivative of neural network wrt its input 
	      const cvm::real factor = input_grad[l+j_elem].item<double>();
	      for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
		  for (size_t l_atom = 0; l_atom < (cv[i_cv]->atom_groups)[k_ag]->size(); ++l_atom) {
		      (*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad = factor_polynomial * factor * (*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad;
		  }
	      }
	  }
      }
  }
}

void colvar::torchANN::apply_force(colvarvalue const &force) {

    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
	// If this CV us explicit gradients, then atomic gradients is already calculated
	// We can apply the force to atom groups directly
	if (cv[i_cv]->is_enabled(f_cvc_explicit_gradient)) {
	    for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
		(cv[i_cv]->atom_groups)[k_ag]->apply_colvar_force(force.real_value);
	    }
	} else {
	    const colvarvalue& current_cv_value = cv[i_cv]->value();
	    colvarvalue cv_force(current_cv_value.type());
	    const cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
	    // get the initial index of this cvc
	    size_t l = cvc_indices[i_cv];
	    for (size_t j_elem = 0; j_elem < current_cv_value.size(); ++j_elem) {
		    cv_force[j_elem] += factor_polynomial * input_grad[l+j_elem].item<double>() * force.real_value;
		}
	    cv[i_cv]->apply_force(cv_force);
	}
    }
}

cvm::real colvar::torchANN::dist2(colvarvalue const &x1, colvarvalue const &x2) const
{
  cvm::real diff = x1.real_value - x2.real_value;
  if (is_enabled(f_cvc_periodic)) 
    diff = (diff < - period * 0.5 ? diff + period : (diff > period * 0.5 ? diff - period : diff));
  return diff * diff;
}

colvarvalue colvar::torchANN::dist2_lgrad(colvarvalue const &x1,
                                          colvarvalue const &x2) const
{
  cvm::real diff = x1.real_value - x2.real_value;
  if (is_enabled(f_cvc_periodic)) 
    diff = (diff < - period * 0.5 ? diff + period : (diff > period * 0.5 ? diff - period : diff));
  return 2.0 * diff;
}


colvarvalue colvar::torchANN::dist2_rgrad(colvarvalue const &x1,
                                          colvarvalue const &x2) const
{
  cvm::real diff = x1.real_value - x2.real_value;
  if (is_enabled(f_cvc_periodic)) 
    diff = (diff < - period * 0.5 ? diff + period : (diff > period * 0.5 ? diff - period : diff));
  return (-2.0) * diff;
}


void colvar::torchANN::wrap(colvarvalue &x_unwrapped) const
{
  if (!is_enabled(f_cvc_periodic)) {
    return;
  }
  cvm::real shift =
    cvm::floor((x_unwrapped.real_value - wrap_center) / period + 0.5);
  x_unwrapped.real_value -= shift * period;
}

#endif
