// src/colvarcomp_harmonicforceconstant.cpp

#include "colvarcomp_harmonicforceconstant.h"
#include "colvarbias_restraint.h" // ���ڿ��������ذ���
#include "colvarmodule.h"

cvc_harmonicforceconstant::cvc_harmonicforceconstant()
    : cvc() // ���������캯������
{
  set_function_type("harmonicForceConstant"); // ���Ҳ��������������

  provide(f_cvc_explicit_gradient, false); // ������ʽ�ݶ�
  provide(f_cvc_gradient, false);          // �����ݶ� (����ʩ����)
  provide(f_cvc_collect_atom_ids, false); // ����ԭ��ID�ռ�

  provide(f_cvc_inv_gradient);             // �ṩ�������� (ͨ�� calc_force_invgrads ʵ��)
  provide(f_cvc_Jacobian);                 // �ṩ�ſɱȵ�������
  
  harmonic_bias = NULL;
  is_linked = false; // ��������ʼ�� is_linked
  k_exponent = 1.0; // Ĭ��ֵ

  // ������'is_extended_lagrangian' ���� CVC �ĳ�Ա��
  // ���� colvar ����ͨ�� f_cv_external ��־���ơ�

  // ������ʹ�û��ຯ�������ñ߽������
  init_scalar_boundaries(0.0, 1.0); 
  
  x.type(colvarvalue::type_scalar);

  // ������'value_' �����ڡ�ʹ�ü̳еĳ�Ա 'x'
  x.real_value = 0.0; // ��ʼ��ֵ (Ĭ�� k=0)

  cvm::log("Initializing a harmonicForceConstant component.\n");
}

int cvc_harmonicforceconstant::init(std::string const &conf)
{
  cvc::init(conf);
  if (!get_keyval(conf, "harmonicName", harmonic_bias_name, std::string(""))) {
    cvm::error("Error: Missing required parameter harmonicName for harmonicForceConstant component.");
    return COLVARS_INPUT_ERROR;
  }
  
  // ����: ���� kExponent
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
  
  // ���� harmonic ƫ�ã���������� CV ����
  // 'parent' ��ӵ����� CVC �� colvar ����������ȷ�ġ�
  harmonic_bias->set_dynamic_k_cv(cv);
  is_linked = true;
  cvm::log("Successfully linked harmonicForceConstant component to harmonic bias '" + harmonic_bias_name + "'.\n");

  return COLVARS_OK;
}

void cvc_harmonicforceconstant::calc_value()
{
  // ֵ�ɻ��������£�����ʲô�������� (��ȷ)
}

void cvc_harmonicforceconstant::calc_gradients()
{
  // ��������������ԭ�����꣬�����ݶ�Ϊ�� (��ȷ)
}

// ������ʵ�� calc_force_invgrads ������ F_lambda
void cvc_harmonicforceconstant::calc_force_invgrads()
{
  if (is_linked && harmonic_bias) {
    // �� bias ��ȡ����һ��ʱ�䲽����� F_lambda
    // ������洢�� CVC �� 'ft' (total force) ��Ա��
    ft.real_value = harmonic_bias->get_k_derivative();
  } else {
    ft.real_value = 0.0;
  }
}

void cvc_harmonicforceconstant::calc_Jacobian_derivative() {
  jd.type(colvarvalue::type_scalar); // ȷ����������
  jd.real_value = 0.0;
}

void cvc_harmonicforceconstant::apply_force(colvarvalue const &force)
{
  // �������˺�������Ϊ�ա�
  // ������ k_cv �ϵ��� (���� ABF, Meta��) �� colvar::update_forces_energy() ����
  // F_lambda (���� harmonic_bias) �� calc_force_invgrads() ���档
  // �� CVC �����κ�ԭ��ֱ�ӽ�����
}