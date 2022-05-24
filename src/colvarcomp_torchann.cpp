#ifdef TORCH

#include "colvarmodule.h"
#include "colvar.h"
#include "colvarvalue.h"
#include "colvarparse.h"
#include "colvarcomp.h"

colvar::torchANN::torchANN(std::string const &conf)
  : cvc(conf)
{
  set_function_type("torchANN");

  if (period != 0.0) {
    enable(f_cvc_periodic);
  }

  if ((wrap_center != 0.0) && !is_enabled(f_cvc_periodic)) {
    cvm::error("Error: wrapAround was defined in a torchANN component,"
                " but its period has not been set.\n");
    return;
  }

  atoms = parse_group(conf, "atoms");

  std::string model_file ;
  get_keyval(conf, "model_file", model_file, std::string(""));
  get_keyval(conf, "m_output_index", m_output_index, 0);

  module = torch::jit::load(model_file);

  torch::Tensor grad ;

  cvm::log("model file name: \"" + model_file + "\".\n" ) ;
}

colvar::torchANN::~torchANN() {
}

void colvar::torchANN::calc_value() {

  std::vector<double> pos_data(atoms->size() * 3) ;
  size_t ia, j;
  // obtain the values of the arguments
  for (ia = 0; ia < atoms->size(); ia++) {
    for (j = 0; j < 3; j++) 
      pos_data[3*ia + j] = (*atoms)[ia].pos[j];
  }
  // change to torch Tensor 
  torch::Tensor arg_tensor = torch::from_blob(pos_data.data(), {1, atoms->size(),3}, torch::TensorOptions().dtype(torch::kFloat64).requires_grad(false));
  std::vector<torch::jit::IValue> inputs = {arg_tensor.to(torch::kFloat32)};

  // evaluate the value of function
  auto outputs = module.forward(inputs).toTensor()[0] ;

  x = outputs[m_output_index].item<double>() ;

  this->wrap(x);
}

void colvar::torchANN::calc_gradients() {

  std::vector<double> pos_data(atoms->size() * 3) ;
  size_t ia, j;

  // obtain the values of the arguments
  for (ia = 0; ia < atoms->size(); ia++) 
    for (j = 0; j < 3; j++) 
      pos_data[3*ia + j] = (*atoms)[ia].pos[j];

  // change to torch Tensor 
  torch::Tensor arg_tensor = torch::from_blob(pos_data.data(), {1, atoms->size(),3}, torch::TensorOptions().dtype(torch::kFloat64).requires_grad(true));
  std::vector<torch::jit::IValue> inputs = {arg_tensor.to(torch::kFloat32)};

  // evaluate the value of function
  auto outputs = module.forward(inputs).toTensor()[0] ;

  outputs[m_output_index].backward({}, false, false);

  torch::Tensor grad = arg_tensor.grad();

  ia = 0 ;

  for (cvm::atom_iter ai = atoms->begin() ; ai != atoms->end(); ai++, ia++) 
  {
    for (size_t j = 0; j < 3; j ++)
      ai->grad[j] = grad[0][ia][j].item<double>() ;
  }
}

void colvar::torchANN::apply_force(colvarvalue const &force) {
  if (!atoms->noforce) {
    atoms->apply_colvar_force(force.real_value);
  }
}

cvm::real colvar::torchANN::dist2(colvarvalue const &x1, colvarvalue const &x2) const
{
  cvm::real diff = x1.real_value - x2.real_value;
  if (is_enabled(f_cvc_periodic)) 
  {
    diff = (diff < - period * 0.5 ? diff + period : (diff > period * 0.5 ? diff - period : diff));
  }
  return diff * diff;
}


colvarvalue colvar::torchANN::dist2_lgrad(colvarvalue const &x1,
                                          colvarvalue const &x2) const
{
  cvm::real diff = x1.real_value - x2.real_value;
  if (is_enabled(f_cvc_periodic)) 
  {
    diff = (diff < - period * 0.5 ? diff + period : (diff > period * 0.5 ? diff - period : diff));
  }
  return 2.0 * diff;
}


colvarvalue colvar::torchANN::dist2_rgrad(colvarvalue const &x1,
                                          colvarvalue const &x2) const
{
  cvm::real diff = x1.real_value - x2.real_value;
  if (is_enabled(f_cvc_periodic)) 
  {
    diff = (diff < - period * 0.5 ? diff + period : (diff > period * 0.5 ? diff - period : diff));
  }
  return (-2.0) * diff;
}


void colvar::torchANN::wrap(colvarvalue &x_unwrapped) const
{
  if ((x_unwrapped.real_value - wrap_center) >= period * 0.5) {
    x_unwrapped.real_value -= period ;
    return;
  }

  if ((x_unwrapped.real_value - wrap_center) < -period * 0.5) {
    x_unwrapped.real_value += period ;
    return;
  }
}

#endif
