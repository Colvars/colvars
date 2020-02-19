#if (__cplusplus >= 201103L)
#include "colvar_neuralnetworkcompute.h"
#include "colvarparse.h"

namespace neuralnetworkCV {
std::map<std::string, std::pair<std::function<double(double)>, std::function<double(double)>>> activation_function_map
{
    {"tanh", {[](double x){return std::tanh(x);}, [](double x){return 1.0 - std::tanh(x) * std::tanh(x);}}},
    {"linear", {[](double x){return x;}, [](double x){return 1.0;}}},
    {"sigmoid", {[](double x){return 1.0 / (1.0 + std::exp(-x));}, [](double x){return std::exp(-x) / ((1.0 + std::exp(-x)) * (1.0 + std::exp(-x)));}}}
};

denseLayer::denseLayer(const std::string& weights_file, const std::string& biases_file, const std::function<double(double)>& f, const std::function<double(double)>& df): m_activation_function(f), m_activation_function_derivative(df) {
    readFromFile(weights_file, biases_file);
}

void denseLayer::readFromFile(const std::string& weights_file, const std::string& biases_file) {
    // parse weights file
    m_weights.clear();
    m_biases.clear();
    std::string line;
    std::ifstream ifs_weights(weights_file.c_str());
    while (std::getline(ifs_weights, line)) {
        std::vector<std::string> splited_data;
        colvarparse::split_string(line, std::string{" "}, splited_data);
        if (splited_data.size() > 0) {
            std::vector<double> weights_tmp(splited_data.size());
            for (size_t i = 0; i < splited_data.size(); ++i) {
                weights_tmp[i] = std::stod(splited_data[i]);
            }
            m_weights.push_back(weights_tmp);
        }
    }
    // parse biases file
    std::ifstream ifs_biases(biases_file.c_str());
    while (std::getline(ifs_biases, line)) {
        std::vector<std::string> splited_data;
        colvarparse::split_string(line, std::string{" "}, splited_data);
        if (splited_data.size() > 0) {
            m_biases.push_back(std::stod(splited_data[0]));
        }
    }
    m_input_size = m_weights[0].size();
    m_output_size = m_weights.size();
}

void denseLayer::setActivationFunction(const std::function<double(double)>& f, const std::function<double(double)>& df) {
    m_activation_function = f;
    m_activation_function_derivative = df;
}

std::vector<double> denseLayer::compute(const std::vector<double>& input) const {
    std::vector<double> output(m_output_size, 0);
    for (size_t i = 0; i < m_output_size; ++i) {
        for (size_t j = 0; j < m_input_size; ++j) {
            output[i] += input[j] * m_weights[i][j];
        }
        output[i] += m_biases[i];
        output[i] = m_activation_function(output[i]);
    }
    return output;
}

std::vector<double> denseLayer::computeTotalDerivative(const std::vector<double>& input) const {
    std::vector<double> output_grad(m_input_size, 0);
    std::vector<double> sum_with_bias(m_output_size, 0);
    for (size_t i = 0; i < m_output_size; ++i) {
        for (size_t j = 0; j < m_input_size; ++j) {
            sum_with_bias[i] += input[j] * m_weights[i][j];
        }
        sum_with_bias[i] += m_biases[i];
    }
    for (size_t j = 0; j < m_input_size; ++j) {
        for (size_t i = 0; i < m_output_size; ++i) {
            output_grad[j] += m_weights[i][j] * m_activation_function_derivative(sum_with_bias[i]);
        }
    }
    return output_grad;
}

double denseLayer::computeGradient(const std::vector<double>& input, const size_t i, const size_t j) const {
    double sum_with_bias = 0;
    for (size_t j_in = 0; j_in < m_input_size; ++j_in) {
        sum_with_bias += input[j_in] * m_weights[i][j_in];
    }
    sum_with_bias += m_biases[i];
    const double grad_ij = m_activation_function_derivative(sum_with_bias) * m_weights[i][j];
    return grad_ij;
}

double denseLayer::computeNumericalGradient(const std::vector<double>& input, const size_t i, const size_t j, const double epsilon) const {
    std::vector<double> input_prev(input);
    std::vector<double> input_next(input);
    input_prev[j] -= epsilon;
    input_next[j] += epsilon;
    std::vector<double> value_prev = compute(input_prev);
    std::vector<double> value_next = compute(input_next);
    const double numerical_gradient = (value_next[i] - value_prev[i]) / (2.0 * epsilon);
    return numerical_gradient;
}

std::vector<std::vector<double>> denseLayer::computeGradient(const std::vector<double>& input) const {
    std::vector<std::vector<double>> output_grad(m_output_size, std::vector<double>(m_input_size, 0));
    for (size_t j = 0; j < m_input_size; ++j) {
        for (size_t i = 0; i < m_output_size; ++i) {
            output_grad[i][j] = computeGradient(input, i, j);
        }
    }
    return output_grad;
}

std::ostream& denseLayer::showInfo(std::ostream& os) {
    os << "Input size: " << m_input_size << std::endl;
    os << "Output size: " << m_output_size << std::endl;
    os << "Weights: \n";
    for (size_t i = 0; i < m_output_size; ++i) {
        for (size_t j = 0; j < m_input_size; ++j) {
            os << m_weights[i][j] << " ";
        }
        os << '\n';
    }
    os << "Biases: \n";
    for (size_t i = 0; i < m_output_size; ++i) {
        os << m_biases[i] << " ";
    }
    os << std::endl;
    return os;
}

neuralNetworkCompute::neuralNetworkCompute(const std::vector<denseLayer>& dense_layers): m_dense_layers(dense_layers) {
    m_layers_output.resize(m_dense_layers.size());
    m_grads_tmp.resize(m_dense_layers.size());
    for (size_t i_layer = 0; i_layer < m_layers_output.size(); ++i_layer) {
        m_layers_output[i_layer].assign(m_dense_layers[i_layer].getOutputSize(), 0);
        m_grads_tmp[i_layer].assign(m_dense_layers[i_layer].getOutputSize(), std::vector<double>(m_dense_layers[i_layer].getInputSize(), 0));
    }
}

bool neuralNetworkCompute::addDenseLayer(const denseLayer& layer) {
    if (m_dense_layers.empty()) {
        // add layer to this ann directly if m_dense_layers is empty
        m_dense_layers.push_back(layer);
        m_layers_output.push_back(std::vector<double>(layer.getOutputSize()));
        m_grads_tmp.push_back(std::vector<std::vector<double>>(layer.getOutputSize(), std::vector<double>(layer.getInputSize(), 0)));
        return true;
    } else {
        // otherwise, we need to check if the output of last layer in m_dense_layers matches the input of layer to be added
        if (m_dense_layers.back().getOutputSize() == layer.getInputSize()) {
            m_dense_layers.push_back(layer);
            m_layers_output.push_back(std::vector<double>(layer.getOutputSize()));
            m_grads_tmp.push_back(std::vector<std::vector<double>>(layer.getOutputSize(), std::vector<double>(layer.getInputSize(), 0)));
            return true;
        } else {
            return false;
        }
    }
}

std::vector<std::vector<double>> neuralNetworkCompute::multiply_matrix(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B) {
    const size_t m = A.size();
    const size_t n = B.size();
    if (A[0].size() != n) {
        std::cerr << "Error on multiplying matrices!\n";
    }
    const size_t t = B[0].size();
    std::vector<std::vector<double>> C(m, std::vector<double>(t, 0.0));
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < t; ++j) {
            for (size_t k = 0; k < n; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

void neuralNetworkCompute::compute() {
    if (m_dense_layers.empty()) {
        return;
    }
    m_layers_output[0] = m_dense_layers[0].compute(m_input);
    for (size_t i_layer = 1; i_layer < m_dense_layers.size(); ++i_layer) {
        m_layers_output[i_layer] = m_dense_layers[i_layer].compute(m_layers_output[i_layer - 1]);
    }
    // gradients of each layer
//     std::vector<std::vector<std::vector<double>>> grads(m_dense_layers.size());
    m_grads_tmp[0] = m_dense_layers[0].computeGradient(m_input);
    for (size_t i_layer = 1; i_layer < m_dense_layers.size(); ++i_layer) {
        m_grads_tmp[i_layer] = m_dense_layers[i_layer].computeGradient(m_layers_output[i_layer - 1]);
    }
    // chain rule
    if (m_dense_layers.size() > 1) {
        m_chained_grad = multiply_matrix(m_grads_tmp[1], m_grads_tmp[0]);
        for (size_t i_layer = 2; i_layer < m_dense_layers.size(); ++i_layer) {
            m_chained_grad = multiply_matrix(m_grads_tmp[i_layer], m_chained_grad);
        }
    } else {
        m_chained_grad = m_grads_tmp[0];
    }
}

std::vector<double> neuralNetworkCompute::compute(const std::vector<double>& input) const {
    if (m_dense_layers.empty()) {
        return input;
    }
    if (input.size() != m_dense_layers.front().getInputSize()) {
        return input;
    }
    // final output
    std::vector<double> output = m_dense_layers[0].compute(input);
    for (size_t i_layer = 1; i_layer < m_dense_layers.size(); ++i_layer) {
        output = m_dense_layers[i_layer].compute(output);
    }
    return output;
}

double neuralNetworkCompute::compute(const std::vector<double>& input, size_t i) const {
    return compute(input)[i];
}

double neuralNetworkCompute::computeGradient(const std::vector<double>& input, const size_t i, const size_t j) const {
    if (m_dense_layers.empty()) {
        return 0.0;
    }
    // input value of each layer
    std::vector<std::vector<double>> values(m_dense_layers.size());
    values[0] = m_dense_layers[0].compute(input);
    for (size_t i_layer = 1; i_layer < m_dense_layers.size(); ++i_layer) {
        values[i_layer] = m_dense_layers[i_layer].compute(values[i_layer - 1]);
    }
    // gradients of each layer
    std::vector<std::vector<std::vector<double>>> grads(m_dense_layers.size());
    grads[0] = m_dense_layers[0].computeGradient(input);
    for (size_t i_layer = 1; i_layer < m_dense_layers.size(); ++i_layer) {
        grads[i_layer] = m_dense_layers[i_layer].computeGradient(values[i_layer - 1]);
    }
    // chain rule
    if (m_dense_layers.size() > 1) {
        std::vector<std::vector<double>> chained_grad = multiply_matrix(grads[1], grads[0]);
        for (size_t i_layer = 2; i_layer < m_dense_layers.size(); ++i_layer) {
            chained_grad = multiply_matrix(grads[i_layer], chained_grad);
        }
        return chained_grad[i][j];
    } else {
        return grads[0][i][j];
    }
}

double neuralNetworkCompute::computeNumericalGradient(const std::vector<double>& input, const size_t i, const size_t j, const double epsilon) const {
    std::vector<double> input_prev(input);
    std::vector<double> input_next(input);
    input_prev[j] -= epsilon;
    input_next[j] += epsilon;
    std::vector<double> value_prev = compute(input_prev);
    std::vector<double> value_next = compute(input_next);
    const double numerical_gradient = (value_next[i] - value_prev[i]) / (2.0 * epsilon);
    return numerical_gradient;
}
}

#endif
