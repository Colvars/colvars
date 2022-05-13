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

  torch::Tensor grad ;

  std::string model_file ;

  get_keyval(conf, "model_file", model_file, std::string(""));

  cvm::log("model file name: \"" + model_file + "\".\n" ) ;
}

colvar::torchANN::~torchANN() {
}

void colvar::torchANN::calc_value() {
}

void colvar::torchANN::calc_gradients() {
}

void colvar::torchANN::apply_force(colvarvalue const &force) {
}

#endif
