// src/colvarcomp_harmonicforceconstant.h

#ifndef COLVARCOMP_HARMONICFORCECONSTANT_H
#define COLVARCOMP_HARMONICFORCECONSTANT_H

#include "colvar.h"
#include "colvarcomp.h"

// ǰ�������Ա���ѭ������
class colvarbias_restraint;

class cvc_harmonicforceconstant : public colvar::cvc {
public:
  cvc_harmonicforceconstant();
  virtual int init(std::string const &conf);
  virtual ~cvc_harmonicforceconstant() {}

  // ��������������� CV �� Bias ��ʼ���󱻵���
  virtual int link_bias(colvarmodule *cvm, colvar *cv);

  virtual void calc_value();
  virtual void calc_gradients();
  
  // ���������� CVC ������ϵͳ����F_lambda������ȷ����
  virtual void calc_force_invgrads(); 

  // �������������Ӧ��Ϊ��
  virtual void apply_force(colvarvalue const &force);

protected:
  // Ҫ���Ƶ� harmonic ƫ�õ�����
  std::string harmonic_bias_name;
  // ָ�� harmonic ƫ��ʵ����ָ��
  colvarbias_restraint *harmonic_bias;
  
  // ��������� is_linked ��Ա����
  bool is_linked;
};

#endif