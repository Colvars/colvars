// src/colvarcomp_harmonicforceconstant.cpp

#include "colvarcomp_harmonicforceconstant.h"
#include "colvarbias_restraint.h" // 现在可以完整地包含
#include "colvarmodule.h"

cvc_harmonicforceconstant::cvc_harmonicforceconstant()
    : cvc() // 修正：构造函数调用
{
  set_function_type("harmonicForceConstant"); // 最好也在这里设置类型

  provide(f_cvc_explicit_gradient, false); // 禁用显式梯度
  provide(f_cvc_gradient, false);          // 禁用梯度 (不能施加力)
  provide(f_cvc_collect_atom_ids, false); // 禁用原子ID收集

  provide(f_cvc_inv_gradient);             // 提供总力计算 (通过 calc_force_invgrads 实现)
  provide(f_cvc_Jacobian);                 // 提供雅可比导数计算
  
  harmonic_bias = NULL;
  is_linked = false; // 修正：初始化 is_linked
  k_exponent = 1.0; // 默认值

  // 修正：'is_extended_lagrangian' 不是 CVC 的成员。
  // 它在 colvar 级别通过 f_cv_external 标志控制。

  // 修正：使用基类函数来设置边界和类型
  init_scalar_boundaries(0.0, 1.0); 
  
  x.type(colvarvalue::type_scalar);

  // 修正：'value_' 不存在。使用继承的成员 'x'
  x.real_value = 0.0; // 初始化值 (默认 k=0)

  cvm::log("Initializing a harmonicForceConstant component.\n");
}

int cvc_harmonicforceconstant::init(std::string const &conf)
{
  cvc::init(conf);
  if (!get_keyval(conf, "harmonicName", harmonic_bias_name, std::string(""))) {
    cvm::error("Error: Missing required parameter harmonicName for harmonicForceConstant component.");
    return COLVARS_INPUT_ERROR;
  }
  
  // 新增: 解析 kExponent
  if (get_keyval(conf, "kExponent", k_exponent, k_exponent)) {
    if (k_exponent <= 0.0) {
        cvm::error("Error: kExponent must be positive for harmonicForceConstant component.", COLVARS_INPUT_ERROR);
        return COLVARS_INPUT_ERROR;
    }
    cvm::log("Using exponent kExponent = " + cvm::to_str(k_exponent) + " for force constant scaling.\n");
  }
  
  return COLVARS_OK;
}

int cvc_harmonicforceconstant::link_bias(colvarmodule *cvm, colvar *cv)
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
  // 'parent' 是拥有这个 CVC 的 colvar 对象，这是正确的。
  harmonic_bias->set_dynamic_k_cv(cv);
  is_linked = true;
  cvm::log("Successfully linked harmonicForceConstant component to harmonic bias '" + harmonic_bias_name + "'.\n");

  return COLVARS_OK;
}

void cvc_harmonicforceconstant::calc_value()
{
  // 值由积分器更新，这里什么都不用做 (正确)
}

void cvc_harmonicforceconstant::calc_gradients()
{
  // 这个组件不依赖于原子坐标，所以梯度为空 (正确)
}

// 修正：实现 calc_force_invgrads 来报告 F_lambda
void cvc_harmonicforceconstant::calc_force_invgrads()
{
  if (is_linked && harmonic_bias) {
    // 从 bias 获取在上一个时间步计算的 F_lambda
    // 并将其存储在 CVC 的 'ft' (total force) 成员中
    ft.real_value = harmonic_bias->get_k_derivative();
  } else {
    ft.real_value = 0.0;
  }
}

void cvc_harmonicforceconstant::calc_Jacobian_derivative() {
  jd.type(colvarvalue::type_scalar); // 确保类型设置
  jd.real_value = 0.0;
}

void cvc_harmonicforceconstant::apply_force(colvarvalue const &force)
{
  // 修正：此函数必须为空。
  // 作用在 k_cv 上的力 (来自 ABF, Meta等) 由 colvar::update_forces_energy() 处理。
  // F_lambda (来自 harmonic_bias) 由 calc_force_invgrads() 报告。
  // 此 CVC 不与任何原子直接交互。
}