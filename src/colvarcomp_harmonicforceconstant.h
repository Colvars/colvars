// src/colvarcomp_harmonicforceconstant.h

#ifndef COLVARCOMP_HARMONICFORCECONSTANT_H
#define COLVARCOMP_HARMONICFORCECONSTANT_H

#include "colvar.h"
#include "colvarcomp.h"

// 前向声明以避免循环依赖
class colvarbias_restraint;

class cvc_harmonicforceconstant : public colvar::cvc {
public:
  cvc_harmonicforceconstant(std::string const &conf);
  virtual int init(std::string const &conf);
  virtual ~cvc_harmonicforceconstant() {}

  // 这个函数将在所有 CV 和 Bias 初始化后被调用
  virtual int link_bias(colvarmodule *cvm);

  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);

protected:
  // 要控制的 harmonic 偏置的名称
  std::string harmonic_bias_name;
  // 指向 harmonic 偏置实例的指针
  colvarbias_restraint *harmonic_bias;
};

#endif
