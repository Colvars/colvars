// src/colvarcomp_harmonicforceconstant.h

#ifndef COLVARCOMP_HARMONICFORCECONSTANT_H
#define COLVARCOMP_HARMONICFORCECONSTANT_H

#include "colvar.h"
#include "colvarcomp.h"

// ǰ�������Ա���ѭ������
class colvarbias_restraint;

class cvc_harmonicforceconstant : public colvar::cvc {
public:
  cvc_harmonicforceconstant(std::string const &conf);
  virtual int init(std::string const &conf);
  virtual ~cvc_harmonicforceconstant() {}

  // ��������������� CV �� Bias ��ʼ���󱻵���
  virtual int link_bias(colvarmodule *cvm);

  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);

protected:
  // Ҫ���Ƶ� harmonic ƫ�õ�����
  std::string harmonic_bias_name;
  // ָ�� harmonic ƫ��ʵ����ָ��
  colvarbias_restraint *harmonic_bias;
};

#endif
