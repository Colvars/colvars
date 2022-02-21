
/// custom expression of colvar components
/// used either in a colvar object or in a customColvar component

class colvar::customFunction
{
protected:
#ifdef LEPTON
  /// Vector of evaluators for custom functions using Lepton
  std::vector<Lepton::CompiledExpression *> value_evaluators;
  /// Vector of evaluators for gradients of custom functions
  std::vector<Lepton::CompiledExpression *> gradient_evaluators;
  /// Vector of references to cvc values to be passed to Lepton evaluators
  std::vector<double *> value_eval_var_refs;
  std::vector<double *> grad_eval_var_refs;
  /// Unused value that is written to when a variable simplifies out of a Lepton expression
  double dev_null;
#endif
public:
  customFunction(std::string const &conf);
  virtual ~customFunction();
  int init_custom_function(std::string const &conf);
  virtual void calc_value();
  virtual void calc_gradients();
  virtual void apply_force(colvarvalue const &force);
};
