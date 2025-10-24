// src/colvarcomp_harmonicforceconstant.cpp

#include "colvarcomp_harmonicforceconstant.h"
#include "colvarbias_restraint.h" // ���ڿ��������ذ���
#include "colvarmodule.h"

cvc_harmonicforceconstant::cvc_harmonicforceconstant(std::string const &conf)
    : cvc(conf)
{
  harmonic_bias = NULL;
  is_extended_lagrangian = true; // ����Ϊ��չ����
  period = 0.0;
  lower_boundary = 0.0; // ֵ������ 0 �� 1 ֮��
  upper_boundary = 1.0;
  value_ = 0.0; // ��ʼ��ֵ
  force_ = 0.0; // ��ʼ����
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
  
  // ���� harmonic ƫ�ã���������� CV ����
  harmonic_bias->set_dynamic_k_cv(parent);
  is_linked = true;
  cvm->log("Successfully linked harmonicForceConstant component to harmonic bias '" + harmonic_bias_name + "'.\n");

  return COLVARS_OK;
}

void cvc_harmonicforceconstant::calc_value()
{
  // ֵ�ɻ��������£�����ʲô��������
}

void cvc_harmonicforceconstant::calc_gradients()
{
  // ��������������ԭ�����꣬�����ݶ�Ϊ��
}

void cvc_harmonicforceconstant::apply_force(colvarvalue const &force)
{
  // �� harmonic ƫ�ý�������ѧ��
  // �������������������������� CV �� "�ٶ�" �� "λ��" (��ֵ)
  force_ += force.real_value;
}
