// src/colvarcomp_harmonicforceconstant.cpp

#include "colvarcomp_harmonicforceconstant.h"
#include "colvarbias_restraint.h" // 现在可以完整地包含
#include "colvarmodule.h"

cvc_harmonicforceconstant::cvc_harmonicforceconstant(std::string const &conf)
    : cvc(conf)
{
  harmonic_bias = NULL;
  is_extended_lagrangian = true; // 声明为扩展坐标
  period = 0.0;
  lower_boundary = 0.0; // 值限制在 0 和 1 之间
  upper_boundary = 1.0;
  value_ = 0.0; // 初始化值
  force_ = 0.0; // 初始化力
  cvm::log("Initializing a harmonicForceConstant component.\n");
}

int cvc_harmonicforceconstant::init(std::string const &conf)
{
  cvc::init(conf);
  if (!get_keyval(conf, "harmonicName", harmonic_bias_name, std::string(""))) {
    cvm->error("Error: Missing required parameter harmonicName for harmonicForceConstant component.");
    return COLVARS_INPUT_ERROR;
  }
  return COLVARS_OK;
}

int cvc_harmonicforceconstant::link_bias(colvarmodule *cvm)
{
  if (is_linked) return COLVARS_OK;

  colvarbias *bias = cvm->bias_by_name(harmonic_bias_name);
  if (!bias) {
    cvm->error("Error: Cannot find harmonic bias named '" + harmonic_bias_name + "' for harmonicForceConstant component.");
    return COLVARS_INPUT_ERROR;
  }
  
  harmonic_bias = dynamic_cast<colvarbias_restraint *>(bias);
  if (!harmonic_bias) {
    cvm->error("Error: Bias '" + harmonic_bias_name + "' is not a harmonic restraint.");
    return COLVARS_INPUT_ERROR;
  }
  
  // 告诉 harmonic 偏置，它将由这个 CV 控制
  harmonic_bias->set_dynamic_k_cv(parent);
  is_linked = true;
  cvm->log("Successfully linked harmonicForceConstant component to harmonic bias '" + harmonic_bias_name + "'.\n");

  return COLVARS_OK;
}

void cvc_harmonicforceconstant::calc_value()
{
  // 值由积分器更新，这里什么都不用做
}

void cvc_harmonicforceconstant::calc_gradients()
{
  // 这个组件不依赖于原子坐标，所以梯度为空
}

void cvc_harmonicforceconstant::apply_force(colvarvalue const &force)
{
  // 从 harmonic 偏置接收热力学力
  // 这个力将被积分器用来更新这个 CV 的 "速度" 和 "位置" (即值)
  force_ += force.real_value;
}
