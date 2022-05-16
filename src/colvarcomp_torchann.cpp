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

  atoms = parse_group(conf, "atoms");

  std::string model_file ;
  get_keyval(conf, "model_file", model_file, std::string(""));

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
  torch::Tensor arg_tensor = torch::from_blob(pos_data.data(), {1, atoms->size(),3}, torch::TensorOptions().dtype(torch::kFloat64).requires_grad(true));
  std::vector<torch::jit::IValue> inputs = {arg_tensor.to(torch::kFloat32)};

  // evaluate the value of function
  auto outputs = module.forward(inputs).toTensor()[0] ;

  x = outputs[0].item<double>() ;
}

void colvar::torchANN::calc_gradients() {

  std::vector<double> pos_data(atoms->size() * 3) ;
  size_t ia, j;
  // obtain the values of the arguments
  for (ia = 0; ia < atoms->size(); ia++) {
    for (j = 0; j < 3; j++) 
      pos_data[3*ia + j] = (*atoms)[ia].pos[j];
  }

  // change to torch Tensor 
  torch::Tensor arg_tensor = torch::from_blob(pos_data.data(), {1, atoms->size(),3}, torch::TensorOptions().dtype(torch::kFloat64).requires_grad(true));
  std::vector<torch::jit::IValue> inputs = {arg_tensor.to(torch::kFloat32)};

  // evaluate the value of function
  auto outputs = module.forward(inputs).toTensor()[0] ;

  outputs[0].backward({}, false, false);

  torch::Tensor grad = arg_tensor.grad();

  ia = 0 ;
  for (cvm::atom_iter ai = atoms->begin() ; ai != atoms->end(); ai++, ia++) 
    for (size_t j = 0; j < 3; j ++)
      ai->grad[j] = grad[0][ia][j].item<double>();
}

void colvar::torchANN::apply_force(colvarvalue const &force) {
  if (!atoms->noforce) {
    atoms->apply_colvar_force(force.real_value);
  }
}

#endif
