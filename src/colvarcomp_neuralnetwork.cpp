#if (__cplusplus >= 201103L)

#include <numeric>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvarparse.h"
#include "colvar.h"
#include "colvarcomp.h"
#include "colvar_neuralnetworkcompute.h"

using namespace neuralnetworkCV;

colvar::customColvar::customColvar(std::string const &conf): linearCombination(conf) {
    use_custom_function = false;
    // code swipe from colvar::init_custom_function
#ifdef LEPTON
    std::string expr_in, expr;
    std::vector<Lepton::ParsedExpression> pexprs;
    Lepton::ParsedExpression pexpr;
    double *ref;
    size_t pos = 0; // current position in config string
    if (key_lookup(conf, "customFunction", &expr_in, &pos)) {
        use_custom_function = true;
        cvm::log("This colvar uses a custom function.\n");
        do {
            expr = expr_in;
            if (cvm::debug())
                cvm::log("Parsing expression \"" + expr + "\".\n");
            try {
                pexpr = Lepton::Parser::parse(expr);
                pexprs.push_back(pexpr);
            } catch (...) {
                cvm::error("Error parsing expression \"" + expr + "\".\n", INPUT_ERROR);
            }
            try {
                value_evaluators.push_back(new Lepton::CompiledExpression(pexpr.createCompiledExpression()));
                // Define variables for cvc values
                for (size_t i = 0; i < cv.size(); ++i) {
                    for (size_t j = 0; j < cv[i]->value().size(); ++j) {
                        std::string vn = cv[i]->name + (cv[i]->value().size() > 1 ? cvm::to_str(j+1) : "");
                        try {
                            ref = &value_evaluators.back()->getVariableReference(vn);
                        } catch (...) {
                            ref = &dev_null;
                            cvm::log("Warning: Variable " + vn + " is absent from expression \"" + expr + "\".\n");
                        }
                        value_eval_var_refs.push_back(ref);
                    }
                }
            } catch (...) {
                cvm::error("Error compiling expression \"" + expr + "\".\n", INPUT_ERROR);
            }
        } while (key_lookup(conf, "customFunction", &expr_in, &pos));
        // Now define derivative with respect to each scalar sub-component
        for (size_t i = 0; i < cv.size(); ++i) {
            for (size_t j = 0; j < cv[i]->value().size(); ++j) {
                std::string vn = cv[i]->name + (cv[i]->value().size() > 1 ? cvm::to_str(j+1) : "");
                for (size_t c = 0; c < pexprs.size(); ++c) {
                    gradient_evaluators.push_back(new Lepton::CompiledExpression(pexprs[c].differentiate(vn).createCompiledExpression()));
                    for (size_t k = 0; k < cv.size(); ++k) {
                        for (size_t l = 0; l < cv[k]->value().size(); l++) {
                            std::string vvn = cv[k]->name + (cv[k]->value().size() > 1 ? cvm::to_str(l+1) : "");
                            try {
                                ref = &gradient_evaluators.back()->getVariableReference(vvn);
                            } catch (...) {
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
        }
        if (value_evaluators.size() != 1) {
            x.type(colvarvalue::type_vector);
        } else {
            x.type(colvarvalue::type_scalar);
        }
    } else {
        cvm::log(std::string{"Warning: no customFunction specified.\n"});
        cvm::log(std::string{"Warning: use linear combination instead.\n"});
    }
#endif
}

colvar::customColvar::~customColvar() {
#ifdef LEPTON
    for (size_t i = 0; i < value_evaluators.size(); ++i) {
        if (value_evaluators[i] != nullptr) delete value_evaluators[i];
    }
    for (size_t i = 0; i < gradient_evaluators.size(); ++i) {
        if (gradient_evaluators[i] != nullptr) delete gradient_evaluators[i];
    }
#endif
}

void colvar::customColvar::calc_value() {
#ifdef LEPTON
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        cv[i_cv]->calc_value();
    }
    x.reset();
    size_t l = 0;
    for (size_t i = 0; i < x.size(); ++i) {
        for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
            colvarvalue current_cv_value(cv[i_cv]->value());
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

void colvar::customColvar::calc_gradients() {
#ifdef LEPTON
    size_t r = 0; // index in the vector of variable references
    size_t e = 0; // index of the gradient evaluator
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) { // for each CV
        cv[i_cv]->calc_gradients();
        if (use_explicit_gradients) {
            colvarvalue current_cv_value(cv[i_cv]->value());
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

void colvar::customColvar::apply_force(colvarvalue const &force) {
#ifdef LEPTON
    size_t r = 0; // index in the vector of variable references
    size_t e = 0; // index of the gradient evaluator
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        // If this CV us explicit gradients, then atomic gradients is already calculated
        // We can apply the force to atom groups directly
        if (use_explicit_gradients) {
            for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
                (cv[i_cv]->atom_groups)[k_ag]->apply_colvar_force(force.real_value);
            }
        } else {
            colvarvalue current_cv_value(cv[i_cv]->value());
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

colvar::neuralNetwork::neuralNetwork(std::string const &conf): linearCombination(conf) {
    function_type = "neuralNetwork";
    // the output of neural network consists of multiple values
    // read "output_component" key to determine it
    get_keyval(conf, "output_component", m_output_index);
    // read weight files
    bool has_weight_files = true;
    size_t num_layers_weight = 0;
    std::vector<std::string> weight_files;
    while (has_weight_files) {
        std::string lookup_key = std::string{"layer"} + cvm::to_str(num_layers_weight + 1) + std::string{"_WeightsFile"};
        if (key_lookup(conf, lookup_key.c_str())) {
            std::string weight_filename;
            get_keyval(conf, lookup_key.c_str(), weight_filename, std::string(""));
            weight_files.push_back(weight_filename);
            cvm::log(std::string{"Will read layer["} + cvm::to_str(num_layers_weight + 1) + std::string{"] weights from "} + weight_filename + '\n');
            ++num_layers_weight;
        } else {
            has_weight_files = false;
        }
    }
    // read bias files
    bool has_bias_files = true;
    size_t num_layers_bias = 0;
    std::vector<std::string> bias_files;
    while (has_bias_files) {
        std::string lookup_key = std::string{"layer"} + cvm::to_str(num_layers_bias + 1) + std::string{"_BiasesFile"};
        if (key_lookup(conf, lookup_key.c_str())) {
            std::string bias_filename;
            get_keyval(conf, lookup_key.c_str(), bias_filename, std::string(""));
            bias_files.push_back(bias_filename);
            cvm::log(std::string{"Will read layer["} + cvm::to_str(num_layers_bias + 1) + std::string{"] biases from "} + bias_filename + '\n');
            ++num_layers_bias;
        } else {
            has_bias_files = false;
        }
    }
    // read activation function strings
    bool has_activation_functions = true;
    size_t num_activation_functions = 0;
    std::vector<std::string> activation_functions;
    while (has_activation_functions) {
        std::string lookup_key = std::string{"layer"} + cvm::to_str(num_activation_functions + 1) + std::string{"_activation"};
        if (key_lookup(conf, lookup_key.c_str())) {
            std::string function_name;
            get_keyval(conf, lookup_key.c_str(), function_name, std::string(""));
            activation_functions.push_back(function_name);
            cvm::log(std::string{"Will read layer["} + cvm::to_str(num_activation_functions + 1) + std::string{"] biases from "} + function_name + '\n');
            ++num_activation_functions;
        } else {
            has_activation_functions = false;
        }
    }
    // expect the three numbers are equal
    if ((num_layers_weight != num_layers_bias) || (num_layers_bias != num_activation_functions)) {
        cvm::error("Error: the numbers of weights, biases and activation functions do not match.\n");
    }
    for (size_t i_layer = 0; i_layer < num_layers_weight; ++i_layer) {
        // query the map of supported activation functions
        const auto& f = activation_function_map[activation_functions[i_layer]].first;
        const auto& df = activation_function_map[activation_functions[i_layer]].second;
        denseLayer d(weight_files[i_layer], bias_files[i_layer], f, df);
        // add a new dense layer to network
        if (nn.addDenseLayer(d)) {
            // show information about the neural network
            cvm::log("Layer " + cvm::to_str(i_layer) + " : has " + cvm::to_str(d.getInputSize()) + " input nodes and " + cvm::to_str(d.getOutputSize()) + " output nodes.\n");
            for (size_t i_output = 0; i_output < d.getOutputSize(); ++i_output) {
                for (size_t j_input = 0; j_input < d.getInputSize(); ++j_input) {
                    cvm::log("    weights[" + cvm::to_str(i_output) + "][" + cvm::to_str(j_input) + "] = " + cvm::to_str(d.getWeight(i_output, j_input)));
                }
                cvm::log("    biases[" + cvm::to_str(i_output) + "] = " + cvm::to_str(d.getBias(i_output)) + "\n");
            }
        } else {
            cvm::error("Error: error on adding a new dense layer.\n");
        }
    }
    nn.input().resize(cv.size());
}

colvar::neuralNetwork::~neuralNetwork() {}

void colvar::neuralNetwork::calc_value() {
    x.reset();
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        cv[i_cv]->calc_value();
        colvarvalue current_cv_value(cv[i_cv]->value());
        // for current nn implementation we have to assume taht types are always scaler
        if (current_cv_value.type() == colvarvalue::type_scalar) {
            nn.input()[i_cv] = cv[i_cv]->sup_coeff * (cvm::pow(current_cv_value.real_value, cv[i_cv]->sup_np));
        } else {
            cvm::error("Error: using of non-scaler component.\n");
        }
    }
    nn.compute();
    x = nn.getOutput(m_output_index);
}

void colvar::neuralNetwork::calc_gradients() {
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        cv[i_cv]->calc_gradients();
        if ( cv[i_cv]->is_enabled(f_cvc_explicit_gradient) &&
            !cv[i_cv]->is_enabled(f_cvc_scalable) &&
            !cv[i_cv]->is_enabled(f_cvc_scalable_com)) {
            const cvm::real factor = nn.getGradient(m_output_index, i_cv);
            const cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
            for (size_t j_elem = 0; j_elem < cv[i_cv]->value().size(); ++j_elem) {
                for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
                    for (size_t l_atom = 0; l_atom < (cv[i_cv]->atom_groups)[k_ag]->size(); ++l_atom) {
                        (*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad = factor_polynomial * factor * (*(cv[i_cv]->atom_groups)[k_ag])[l_atom].grad;
                    }
                }
            }
        }
    }
}

void colvar::neuralNetwork::apply_force(colvarvalue const &force) {
    for (size_t i_cv = 0; i_cv < cv.size(); ++i_cv) {
        // If this CV us explicit gradients, then atomic gradients is already calculated
        // We can apply the force to atom groups directly
        if ( cv[i_cv]->is_enabled(f_cvc_explicit_gradient) &&
            !cv[i_cv]->is_enabled(f_cvc_scalable) &&
            !cv[i_cv]->is_enabled(f_cvc_scalable_com)
        ) {
            for (size_t k_ag = 0 ; k_ag < cv[i_cv]->atom_groups.size(); ++k_ag) {
                (cv[i_cv]->atom_groups)[k_ag]->apply_colvar_force(force.real_value);
            }
        } else {
            // Compute factors for polynomial combinations
            const cvm::real factor_polynomial = getPolynomialFactorOfCVGradient(i_cv);
            const cvm::real factor = nn.getGradient(m_output_index, i_cv);;
            colvarvalue cv_force = force.real_value * factor * factor_polynomial;
            cv[i_cv]->apply_force(cv_force);
        }
    }
}

#endif
