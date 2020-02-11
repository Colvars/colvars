#if (__cplusplus >= 201103L)
#ifndef NEURALNETWORKCOMPUTE_H
#define NEURALNETWORKCOMPUTE_H

#include <vector>
#include <functional>
#include <string>
#include <iostream>
#include <cmath>
#include <map>

namespace neuralnetworkCV {
/// mapping from a string to the activation function and its derivative
extern std::map<std::string, std::pair<std::function<double(double)>, std::function<double(double)>>> activation_function_map;

class denseLayer {
private:
    size_t m_input_size;
    size_t m_output_size;
    std::function<double(double)> m_activation_function;
    std::function<double(double)> m_activation_function_derivative;
    /// weights[i][j] is the weight of the i-th output and the j-th input
    std::vector<std::vector<double>> m_weights;
    /// bias of each node
    std::vector<double> m_biases;
public:
    /// empty constructor
    denseLayer() {}
    /*! @param[in]  weights_file    filename of the weights file
     *  @param[in]  biases_file     filename of the biases file
     *  @param[in]  f               activation function
     *  @param[in]  df              derivative of the activation function
     */
    denseLayer(const std::string& weights_file, const std::string& biases_file, const std::function<double(double)>& f, const std::function<double(double)>& df);
    /// read data from file
    void readFromFile(const std::string& weights_file, const std::string& biases_file);
    /// setup activation function
    void setActivationFunction(const std::function<double(double)>& f, const std::function<double(double)>& df);
    /// compute the value of this layer
    std::vector<double> compute(const std::vector<double>& input) const;
    /// compute the total derivative of this layer wrt input
    std::vector<double> computeTotalDerivative(const std::vector<double>& input) const;
    /// compute the gradient of i-th output wrt j-th input
    double computeGradient(const std::vector<double>& input, const size_t i, const size_t j) const;
    /// compute the numerical gradient of i-th output wrt j-th input
    double computeNumericalGradient(const std::vector<double>& input, const size_t i, const size_t j, const double epsilon = 0.0001) const;
    /// output[i][j] is the gradient of i-th output wrt j-th input
    std::vector<std::vector<double>> computeGradient(const std::vector<double>& input) const;
    /// dump info
    std::ostream& showInfo(std::ostream& os = std::cout);
    /// get the input size
    size_t getInputSize() const {
        return m_input_size;
    }
    /// get the output size
    size_t getOutputSize() const {
        return m_output_size;
    }
    /// getter for weights and biases
    double getWeight(size_t i, size_t j) const {
        return m_weights[i][j];
    }
    double getBias(size_t i) const {
        return m_biases[i];
    }
    ~denseLayer() {}
};

class neuralNetworkCompute {
private:
    std::vector<denseLayer> m_dense_layers;
private:
    /// helper function: multiply two matrix constructed from 2D vector
    static std::vector<std::vector<double>> multiply_matrix(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B);
public:
    neuralNetworkCompute(): m_dense_layers(0) {}
    neuralNetworkCompute(const std::vector<denseLayer>& dense_layers): m_dense_layers(dense_layers) {}
    bool addDenseLayer(const denseLayer& layer);
    /// compute the values of all output nodes
    std::vector<double> compute(const std::vector<double>& input) const;
    /// compute the value of a specified output node
    double compute(const std::vector<double>& input, size_t i) const;
    /// compute the gradient of i-th output wrt j-th input
    double computeGradient(const std::vector<double>& input, const size_t i, const size_t j) const;
    /// compute the numerical gradient of i-th output wrt j-th input
    double computeNumericalGradient(const std::vector<double>& input, const size_t i, const size_t j, const double epsilon = 0.0001) const;
    /// get a specified layer
    const denseLayer& getLayer(const size_t i) const {return m_dense_layers[i];}
    /// get the number of layers
    size_t getNumberOfLayers() const {return m_dense_layers.size();}
};

}
#endif
#endif
