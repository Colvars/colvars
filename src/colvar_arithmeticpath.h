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
T logsumexp(const vector<T>& a, const vector<T>& b, T* sign_factor = nullptr) {
    const auto max_a = *std::max_element(a.begin(), a.end());
    T sum = T();
    for (size_t i = 0; i < a.size(); ++i) {
        sum += b[i] * cvm::exp(a[i] - max_a);
    }
    if (sign_factor) {
        if (sum > 0) {
            *sign_factor = 1.0;
        } else if (sum < 0) {
            sum = -sum;
            *sign_factor = -1.0;
        } else {
            sum = 1.0;
            *sign_factor = 0.0;
        }
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


enum path_sz {S, Z};

template <typename element_type, typename scalar_type, path_sz path_type>
class ArithmeticPathBase {
public:
    ArithmeticPathBase() {}
    virtual ~ArithmeticPathBase() {}
    virtual void initialize(size_t p_num_elements, size_t p_total_frames, double p_lambda, const vector<element_type>& p_element, const vector<double>& p_weights);
    virtual void updateDistanceToReferenceFrames() = 0;
    virtual void computeValue();
    virtual void computeDerivatives();
    virtual void compute();
    virtual void reComputeLambda(const vector<scalar_type>& rmsd_between_refs);
protected:
    scalar_type lambda;
    vector<scalar_type> weights;
    vector<scalar_type> frame_indexes;
    size_t num_elements;
    size_t total_frames;
    vector< vector<element_type> > frame_element_distances;
    scalar_type s;
    scalar_type z;
    vector<element_type> dsdx;
    vector<element_type> dzdx;
private:
    // intermediate variables
    vector<scalar_type> exponents;
    vector<scalar_type> exponents2;
    vector<scalar_type> exponents3;
    scalar_type normalization_factor;
    scalar_type log_sum_exp_0;
    scalar_type log_sum_exp_1;
};

template <typename element_type, typename scalar_type, path_sz path_type>
void ArithmeticPathBase<element_type, scalar_type, path_type>::initialize(size_t p_num_elements, size_t p_total_frames, double p_lambda, const vector<element_type>& p_element, const vector<double>& p_weights) {
    lambda = p_lambda;
    weights = p_weights;
    num_elements = p_num_elements;
    total_frames = p_total_frames;
    frame_element_distances.resize(total_frames, p_element);
    frame_indexes.resize(total_frames);
    exponents.resize(total_frames);
    exponents2.resize(total_frames);
    exponents3.resize(total_frames);
    for (size_t i_frame = 0; i_frame < frame_element_distances.size(); ++i_frame) {
        for (size_t j_elem = 0; j_elem < num_elements; ++j_elem) {
            frame_element_distances[i_frame][j_elem].reset();
        }
        frame_indexes[i_frame] = i_frame;
    }
    s = scalar_type(0);
    z = scalar_type(0);
    dsdx = p_element;
    dzdx = p_element;
    log_sum_exp_0 = scalar_type(0);
    log_sum_exp_1 = scalar_type(0);
    normalization_factor = 1.0 / static_cast<scalar_type>(total_frames - 1);
}

template <typename element_type, typename scalar_type, path_sz path_type>
void ArithmeticPathBase<element_type, scalar_type, path_type>::computeValue() {
    updateDistanceToReferenceFrames();
    for (size_t i_frame = 0; i_frame < frame_element_distances.size(); ++i_frame) {
        scalar_type exponent_tmp = scalar_type(0);
        for (size_t j_elem = 0; j_elem < num_elements; ++j_elem) {
            exponent_tmp += weights[j_elem] * frame_element_distances[i_frame][j_elem] * weights[j_elem] * frame_element_distances[i_frame][j_elem];
        }
        exponents[i_frame] = exponent_tmp * -1.0 * lambda;
    }
    log_sum_exp_0 = logsumexp(exponents);
    log_sum_exp_1 = logsumexp(exponents, frame_indexes);
    s = normalization_factor * cvm::exp(log_sum_exp_1 - log_sum_exp_0);
    z = -1.0 / lambda * log_sum_exp_0;
}

template <typename element_type, typename scalar_type, path_sz path_type>
void ArithmeticPathBase<element_type, scalar_type, path_type>::compute() {
    computeValue();
    computeDerivatives();
}

template <typename element_type, typename scalar_type, path_sz path_type>
void ArithmeticPathBase<element_type, scalar_type, path_type>::computeDerivatives() {
    for (size_t j_elem = 0; j_elem < num_elements; ++j_elem) {
        for (size_t i_frame = 0; i_frame < frame_element_distances.size(); ++i_frame) {
            exponents2[i_frame] = frame_element_distances[i_frame][j_elem];
            exponents3[i_frame] = static_cast<scalar_type>(i_frame) * frame_element_distances[i_frame][j_elem];
        }
        scalar_type sign_factor_2, sign_factor_3;
        const auto log_sum_exp_2 = logsumexp(exponents, exponents2, &sign_factor_2);
        const auto log_sum_exp_3 = logsumexp(exponents, exponents3, &sign_factor_3);
        const auto tmp_factor = 2.0 * weights[j_elem] * weights[j_elem];
        dsdx[j_elem] = -tmp_factor * lambda * normalization_factor *
                        (sign_factor_3 * cvm::exp(log_sum_exp_3 - log_sum_exp_0) - sign_factor_2 * cvm::exp(log_sum_exp_2 + log_sum_exp_1 - 2.0 * log_sum_exp_0));
        dzdx[j_elem] = tmp_factor * sign_factor_2 * cvm::exp(log_sum_exp_2 - log_sum_exp_0);
    }
}

template <typename element_type, typename scalar_type, path_sz path_type>
void ArithmeticPathBase<element_type, scalar_type, path_type>::reComputeLambda(const vector<scalar_type>& rmsd_between_refs) {
    scalar_type mean_square_displacements = 0.0;
    for (size_t i_frame = 1; i_frame < total_frames; ++i_frame) {
        cvm::log(std::string("Distance between frame ") + cvm::to_str(i_frame) + " and " + cvm::to_str(i_frame + 1) + " is " + cvm::to_str(rmsd_between_refs[i_frame - 1]) + std::string("\n"));
        mean_square_displacements += rmsd_between_refs[i_frame - 1] * rmsd_between_refs[i_frame - 1];
    }
    mean_square_displacements /= scalar_type(total_frames - 1);
    lambda = 1.0 / mean_square_displacements;
}
}

#endif // ARITHMETICPATHCV_H
