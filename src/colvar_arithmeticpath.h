#ifndef ARITHMETICPATHCV_H
#define ARITHMETICPATHCV_H

#include "colvarmodule.h"

#include <vector>
#include <cmath>
#include <limits>
#include <string>
#include <algorithm>

namespace ArithmeticPathCV {

using std::vector;

template <typename T>
T logsumexp(const vector<T>& a, const vector<T>& b) {
    const auto max_a = *std::max_element(a.begin(), a.end());
    T sum = T();
    for (size_t i = 0; i < a.size(); ++i) {
        sum += b[i] * cvm::exp(a[i] - max_a);
    }
    return max_a + cvm::logn(sum);
}

template <typename T>
T logsumexp(const vector<T>& a) {
    const auto max_a = *std::max_element(a.begin(), a.end());
    T sum = T();
    for (size_t i = 0; i < a.size(); ++i) {
        sum += cvm::exp(a[i] - max_a);
    }
    return max_a + cvm::logn(sum);
}


template <typename scalar_type>
class ArithmeticPathBase {
public:
    ArithmeticPathBase() {}
    ~ArithmeticPathBase() {}
    void initialize(size_t p_num_elements, size_t p_total_frames, scalar_type p_lambda, const vector<scalar_type>& p_weights);
    void reComputeLambda(const vector<scalar_type>& rmsd_between_refs);
    template <typename element_type>
    void computeValue(const vector<vector<element_type>>& frame_element_distances, scalar_type *s = nullptr, scalar_type *z = nullptr);
    // can only be called after computeValue() for element-wise derivatives and store derivatives of i-th frame to dsdx and dzdx
    template <typename element_type>
    void computeDerivatives(const vector<vector<element_type>>& frame_element_distances, vector<vector<element_type>> *dsdx = nullptr, vector<vector<element_type>> *dzdx = nullptr);
protected:
    scalar_type lambda;
    vector<scalar_type> squared_weights;
    vector<scalar_type> frame_indexes;
    size_t num_elements;
    size_t total_frames;
    vector<scalar_type> exponents;
    scalar_type normalization_factor;
    scalar_type saved_s;
};

template <typename scalar_type>
void ArithmeticPathBase<scalar_type>::initialize(size_t p_num_elements, size_t p_total_frames, scalar_type p_lambda, const vector<scalar_type>& p_weights) {
    lambda = p_lambda;
    for (size_t i = 0; i < p_weights.size(); ++i) squared_weights.push_back(p_weights[i] * p_weights[i]);
    num_elements = p_num_elements;
    total_frames = p_total_frames;
    frame_indexes.resize(total_frames);
    for (size_t i = 0; i < total_frames; ++i) frame_indexes[i] = i;
    exponents.resize(total_frames);
    normalization_factor = 1.0 / static_cast<scalar_type>(total_frames - 1);
    saved_s = scalar_type();
}

template <typename scalar_type>
template <typename element_type>
void ArithmeticPathBase<scalar_type>::computeValue(
    const vector<vector<element_type>>& frame_element_distances,
    scalar_type *s, scalar_type *z)
{
    for (size_t i_frame = 0; i_frame < total_frames; ++i_frame) {
        scalar_type exponent_tmp = scalar_type();
        for (size_t j_elem = 0; j_elem < num_elements; ++j_elem) {
            exponent_tmp += squared_weights[j_elem] * frame_element_distances[i_frame][j_elem] * frame_element_distances[i_frame][j_elem];
        }
        exponents[i_frame] = exponent_tmp * -1.0 * lambda;
    }
    const scalar_type log_sum_exp_0 = logsumexp(exponents);
    const scalar_type log_sum_exp_1 = logsumexp(exponents, frame_indexes);
    saved_s = normalization_factor * cvm::exp(log_sum_exp_1 - log_sum_exp_0);
    if (s != nullptr) {
        *s = saved_s;
    }
    if (z != nullptr) {
        *z = -1.0 / lambda * log_sum_exp_0;
    }
}

template <typename scalar_type>
void ArithmeticPathBase<scalar_type>::reComputeLambda(const vector<scalar_type>& rmsd_between_refs) {
    scalar_type mean_square_displacements = 0.0;
    for (size_t i_frame = 1; i_frame < total_frames; ++i_frame) {
        cvm::log(std::string("Distance between frame ") + cvm::to_str(i_frame) + " and " + cvm::to_str(i_frame + 1) + " is " + cvm::to_str(rmsd_between_refs[i_frame - 1]) + std::string("\n"));
        mean_square_displacements += rmsd_between_refs[i_frame - 1] * rmsd_between_refs[i_frame - 1];
    }
    mean_square_displacements /= scalar_type(total_frames - 1);
    lambda = 1.0 / mean_square_displacements;
}

// frame-wise derivatives for frames using optimal rotation
template <typename scalar_type>
template <typename element_type>
void ArithmeticPathBase<scalar_type>::computeDerivatives(
    const vector<vector<element_type>>& frame_element_distances,
    vector<vector<element_type>> *dsdx,
    vector<vector<element_type>> *dzdx)
{
    const vector<scalar_type> softmax_out = softmax(exponents);
    vector<scalar_type> tmps;
    tmps.reserve(total_frames);
    for (size_t i_frame = 0; i_frame < total_frames; ++i_frame) {
        tmps.push_back(
            (static_cast<scalar_type>(i_frame) -
             static_cast<scalar_type>(total_frames - 1) * saved_s) *
             normalization_factor);
    }
    if (dsdx != nullptr) {
        for (size_t i_frame = 0; i_frame < total_frames; ++i_frame) {
            for (size_t j_elem = 0; j_elem < num_elements; ++j_elem) {
                (*dsdx)[i_frame][j_elem] =
                    -2.0 * squared_weights[j_elem] * lambda *
                    frame_element_distances[i_frame][j_elem] *
                    softmax_out[i_frame] * tmps[i_frame];
            }
        }
    }
    if (dzdx != nullptr) {
        for (size_t i_frame = 0; i_frame < total_frames; ++i_frame) {
            for (size_t j_elem = 0; j_elem < num_elements; ++j_elem) {
                (*dzdx)[i_frame][j_elem] =
                    2.0 * squared_weights[j_elem] * softmax_out[i_frame] *
                    frame_element_distances[i_frame][j_elem];
            }
        }
    }
}
}

#endif // ARITHMETICPATHCV_H
