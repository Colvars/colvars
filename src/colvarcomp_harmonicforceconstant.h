// src/colvarcomp_harmonicforceconstant.h

#ifndef COLVARCOMP_HARMONICFORCECONSTANT_H
#define COLVARCOMP_HARMONICFORCECONSTANT_H

#include "colvar.h"
#include "colvarcomp.h"

// 前向声明以避免循环依赖
class colvarbias_restraint;

class cvc_harmonicforceconstant : public colvar::cvc {
public:
  cvc_harmonicforceconstant();
  virtual int init(std::string const &conf);
  virtual ~cvc_harmonicforceconstant() {}

  // 这个函数将在所有 CV 和 Bias 初始化后被调用
  virtual int link_bias(colvarmodule *cvm, colvar *cv);

  virtual void calc_value();
  virtual void calc_gradients();
  
  // 修正：这是 CVC 报告其系统力（F_lambda）的正确函数
  virtual void calc_force_invgrads(); 

  // 修正：这个函数应该为空
  virtual void apply_force(colvarvalue const &force);

protected:
  // 要控制的 harmonic 偏置的名称
  std::string harmonic_bias_name;
  // 指向 harmonic 偏置实例的指针
  colvarbias_restraint *harmonic_bias;
  
  // 修正：添加 is_linked 成员变量
  bool is_linked;
};

#endif