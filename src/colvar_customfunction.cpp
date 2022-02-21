
colvar::customFunction::customFunction(std::string const &conf)
{
  init_custom_function(conf);
}



colvar::customFunction::~customFunction() {
#ifdef LEPTON
  for (std::vector<Lepton::CompiledExpression *>::iterator cei = value_evaluators.begin();
       cei != value_evaluators.end();
       ++cei) {
    if (*cei != NULL) delete (*cei);
  }
  value_evaluators.clear();

  for (std::vector<Lepton::CompiledExpression *>::iterator gei = gradient_evaluators.begin();
       gei != gradient_evaluators.end();
       ++gei) {
    if (*gei != NULL) delete (*gei);
  }
  gradient_evaluators.clear();
#endif
}



#ifdef LEPTON
int colvar::init_custom_function(std::string const &conf)
{
  std::string expr, expr_in; // expr_in is a buffer to remember expr after unsuccessful parsing
  std::vector<Lepton::ParsedExpression> pexprs;
  Lepton::ParsedExpression pexpr;
  size_t pos = 0; // current position in config string
  double *ref;

  if (!key_lookup(conf, "customFunction", &expr_in, &pos)) {
    return COLVARS_OK;
  }

  cvm::main()->cite_feature("Custom functions (Lepton)");

  enable(f_cv_custom_function);
  cvm::log("This colvar uses a custom function.\n");

  do {
    expr = expr_in;
    if (cvm::debug())
      cvm::log("Parsing expression \"" + expr + "\".\n");
    try {
      pexpr = Lepton::Parser::parse(expr);
      pexprs.push_back(pexpr);
    }
    catch (...) {
      cvm::error("Error parsing expression \"" + expr + "\".\n", INPUT_ERROR);
      return INPUT_ERROR;
    }

    try {
      value_evaluators.push_back(
          new Lepton::CompiledExpression(pexpr.createCompiledExpression()));
      // Define variables for cvc values
      // Stored in order: expr1, cvc1, cvc2, expr2, cvc1...
      for (size_t i = 0; i < cvcs.size(); i++) {
        for (size_t j = 0; j < cvcs[i]->value().size(); j++) {
          std::string vn = cvcs[i]->name +
              (cvcs[i]->value().size() > 1 ? cvm::to_str(j+1) : "");
          try {
            ref =&value_evaluators.back()->getVariableReference(vn);
          }
          catch (...) { // Variable is absent from expression
            // To keep the same workflow, we use a pointer to a double here
            // that will receive CVC values - even though none was allocated by Lepton
            ref = &dev_null;
            cvm::log("Warning: Variable " + vn + " is absent from expression \"" + expr + "\".\n");
          }
          value_eval_var_refs.push_back(ref);
        }
      }
    }
    catch (...) {
      cvm::error("Error compiling expression \"" + expr + "\".\n", INPUT_ERROR);
      return INPUT_ERROR;
    }
  } while (key_lookup(conf, "customFunction", &expr_in, &pos));


  // Now define derivative with respect to each scalar sub-component
  for (size_t i = 0; i < cvcs.size(); i++) {
    for (size_t j = 0; j < cvcs[i]->value().size(); j++) {
      std::string vn = cvcs[i]->name +
          (cvcs[i]->value().size() > 1 ? cvm::to_str(j+1) : "");
      // Element ordering: we want the
      // gradient vector of derivatives of all elements of the colvar
      // wrt to a given element of a cvc ([i][j])
      for (size_t c = 0; c < pexprs.size(); c++) {
        gradient_evaluators.push_back(
            new Lepton::CompiledExpression(pexprs[c].differentiate(vn).createCompiledExpression()));
        // and record the refs to each variable in those expressions
        for (size_t k = 0; k < cvcs.size(); k++) {
          for (size_t l = 0; l < cvcs[k]->value().size(); l++) {
            std::string vvn = cvcs[k]->name +
                (cvcs[k]->value().size() > 1 ? cvm::to_str(l+1) : "");
            try {
              ref = &gradient_evaluators.back()->getVariableReference(vvn);
            }
            catch (...) { // Variable is absent from derivative
              // To keep the same workflow, we use a pointer to a double here
              // that will receive CVC values - even though none was allocated by Lepton
              cvm::log("Warning: Variable " + vvn + " is absent from derivative of \"" + expr + "\" wrt " + vn + ".\n");
              ref = &dev_null;
            }
            grad_eval_var_refs.push_back(ref);
          }
        }
      }
    }
  }


  if (value_evaluators.size() == 0) {
    cvm::error("Error: no custom function defined.\n", INPUT_ERROR);
    return INPUT_ERROR;
  }

  std::string type_str;
  bool b_type_specified = get_keyval(conf, "customFunctionType",
                                     type_str, "scalar", parse_silent);
  x.type(colvarvalue::type_notset);
  int t;
  for (t = 0; t < colvarvalue::type_all; t++) {
    if (type_str == colvarvalue::type_keyword(colvarvalue::Type(t))) {
      x.type(colvarvalue::Type(t));
      break;
    }
  }
  if (x.type() == colvarvalue::type_notset) {
    cvm::error("Could not parse custom colvar type.", INPUT_ERROR);
    return INPUT_ERROR;
  }

  // Guess type based on number of expressions
  if (!b_type_specified) {
    if (value_evaluators.size() == 1) {
      x.type(colvarvalue::type_scalar);
    } else {
      x.type(colvarvalue::type_vector);
    }
  }

  if (x.type() == colvarvalue::type_vector) {
    x.vector1d_value.resize(value_evaluators.size());
  }

  x_reported.type(x);
  cvm::log(std::string("Expecting colvar value of type ")
    + colvarvalue::type_desc(x.type())
    + (x.type()==colvarvalue::type_vector ? " of size " + cvm::to_str(x.size()) : "")
    + ".\n");

  if (x.size() != value_evaluators.size()) {
    cvm::error("Error: based on custom function type, expected "
               + cvm::to_str(x.size()) + " scalar expressions, but "
               + cvm::to_str(value_evaluators.size()) + " were found.\n");
    return INPUT_ERROR;
  }

  return COLVARS_OK;
}

#else

int colvar::init_custom_function(std::string const &conf)
{

  std::string expr;
  size_t pos = 0;
  if (key_lookup(conf, "customFunction", &expr, &pos)) {
    std::string msg("Error: customFunction requires the Lepton library.");
#if (__cplusplus < 201103L)
    // NOTE: this is not ideal; testing for the Lepton library's version would
    // be more accurate, but also less portable
    msg +=
      std::string("  Note also that recent versions of Lepton require C++11: "
                  "please see https://colvars.github.io/README-c++11.html.");
#endif
    return cvm::error(msg, COLVARS_NOT_IMPLEMENTED);
  }

  return COLVARS_OK;
}

#endif // #ifdef LEPTON




void colvar::customFunction::calc_value() {
#ifdef LEPTON
  for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
    cv[i_cv]->calc_value();
  }
  x.reset();
  size_t l = 0;
  for (size_t i = 0; i < x.size(); ++i) {
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
      const colvarvalue& current_cv_value = cv[i_cv]->value();
      for (size_t j_elem = 0; j_elem < current_cv_value.size(); ++j_elem) {
        if (current_cv_value.type() == colvarvalue::type_scalar) {
          *(value_eval_var_refs[l++]) = cv[i_cv]->sup_coeff * (cvm::pow(current_cv_value.real_value, cv[i_cv]->sup_np));
        } else {
          *(value_eval_var_refs[l++]) = cv[i_cv]->sup_coeff * current_cv_value[j_elem];
        }
      }
    }
    x[i] = value_evaluators[i]->evaluate();
  }
#endif
  if (!use_custom_function) {
    colvar::linearCombination::calc_value();
  }
}

void colvar::customFunction::calc_gradients() {
#ifdef LEPTON
  size_t r = 0; // index in the vector of variable references
  size_t e = 0; // index of the gradient evaluator
  for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) { // for each CV
    cv[i_cv]->calc_gradients();
    if (cv[i_cv]->is_enabled(f_cvc_explicit_gradient)) {
      const colvarvalue& current_cv_value = cv[i_cv]->value();
      const cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
      for (size_t j_elem = 0; j_elem < current_cv_value.size(); ++j_elem) { // for each element in this CV
        for (size_t c = 0; c < x.size(); ++c) { // for each custom function expression
          for (size_t k = 0; k < cv.size(); ++k) { // this is required since we need to feed all CV values to this expression
            const cvm::real factor_polynomial_k = getPolynomialFactorOfCVGradient(k);
            for (size_t l = 0; l < cv[k]->value().size(); ++l) {
              *(grad_eval_var_refs[r++]) = factor_polynomial_k * cv[k]->value()[l];
            }
          }
          const double expr_grad = gradient_evaluators[e++]->evaluate();
          for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
            for (size_t l_atom = 0; l_atom < (cv[i_cv]->atom_groups)[k_ag]->size(); ++l_atom) {
              (*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad = expr_grad * factor_polynomial * (*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad;
            }
          }
        }
      }
    }
  }
#endif
  if (!use_custom_function) {
    colvar::linearCombination::calc_gradients();
  }
}

void colvar::customFunction::apply_force(colvarvalue const &force) {
#ifdef LEPTON
  size_t r = 0; // index in the vector of variable references
  size_t e = 0; // index of the gradient evaluator
  for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
    // If this CV us explicit gradients, then atomic gradients is already calculated
    // We can apply the force to atom groups directly
    if (cv[i_cv]->is_enabled(f_cvc_explicit_gradient)) {
      for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
        (cv[i_cv]->atom_groups)[k_ag]->apply_colvar_force(force.real_value);
      }
    } else {
      const colvarvalue& current_cv_value = cv[i_cv]->value();
      colvarvalue cv_force(current_cv_value.type());
      const cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
      for (size_t j_elem = 0; j_elem < current_cv_value.size(); ++j_elem) {
        for (size_t c = 0; c < x.size(); ++c) {
          for (size_t k = 0; k < cv.size(); ++k) {
            const cvm::real factor_polynomial_k = getPolynomialFactorOfCVGradient(k);
            for (size_t l = 0; l < cv[k]->value().size(); ++l) {
              *(grad_eval_var_refs[r++]) = factor_polynomial_k * cv[k]->value()[l];
            }
          }
          cv_force[j_elem] += factor_polynomial * gradient_evaluators[e++]->evaluate() * force.real_value;
        }
      }
      cv[i_cv]->apply_force(cv_force);
    }
  }
#endif
  if (!use_custom_function) {
    colvar::linearCombination::apply_force(force);
  }
}

#endif // __cplusplus >= 201103L