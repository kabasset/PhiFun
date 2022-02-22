// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef _PHIMILLER_LAMBDAINTEGRATOR_H
#define _PHIMILLER_LAMBDAINTEGRATOR_H

#include <complex>
#include <vector>

namespace Phi {
namespace Miller {

/**
 * @brief A simplistic tridiagonal matrix.
 * @details
 * Could be vastly optimized, and Eigen probably provides something.
 */
struct TridiagonalMatrix {
  TridiagonalMatrix(std::size_t size) : upper(size - 1), diagonal(size), lower(size - 1) {}
  std::vector<double> upper;
  std::vector<double> diagonal;
  std::vector<double> lower;
};

/**
 * @brief Spline interpolation and integration as proposed by Lance Miller.
 * @details
 * The monochromatic TFs are interpolated along lambdas at a few nanometers,
 * and then integrated along lambdas.
 * Lance Miller proposes to perform both operations at once,
 * i.e. for each position, compute the weight of interpolating splines and perform the summation,
 * 
 * For example, go directly from 40 monochromatic TFs to 1 broadband TF
 * instead of computing 200 interpolated TFs in-between.
 * 
 * @see https://gitlab.euclid-sgs.uk/PF-SHE/SHE_PSFToolkit/-/blob/develop/SHE_PSFToolkit/python/SHE_PSFToolkit/broadband/basis_interpolation.py#L364
 * @see https://math.stackexchange.com/questions/62360/natural-cubic-splines-vs-piecewise-hermite-splines/62454#62454
 * @see https://math.stackexchange.com/questions/3135857/not-a-knot-cubic-spline-interpolation-using-tridiagonal-solver
 * @see https://www.math.uh.edu/~jingqiu/math4364/spline.pdf
 * @see https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
 */
class SplineIntegrator {

public:
  /**
   * @brief Constructor.
   * @param u The input knots x's
   * @param x The output x's
   * @details
   * Precompute spline coefficients.
   */
  SplineIntegrator(std::vector<double> u, std::vector<double> x) :
      m_size(u.size() - 1), m_u(std::move(u)), m_h(m_size), m_x(std::move(x)), m_i(m_size) {
    std::size_t j = 0;
    for (std::size_t i = 0; i < m_size; ++i) {
      m_h[i] = m_u[i + 1] - m_u[i];
      while (m_x[j] < m_u[i + 1]) {
        ++m_i[i];
        const double xjui0 = m_x[j] - m_u[i];
        const double ui1xj = m_u[i + 1] - m_x[j];
        m_k[j].y0 = ui1xj / m_h[i];
        m_k[j].y1 = xjui0 / m_h[i];
        m_k[j].z0 = ui1xj * ui1xj * ui1xj / (6. * m_h[i]) - m_h[i] * ui1xj / 6.;
        m_k[j].z1 = xjui0 * xjui0 * xjui0 / m_h[i] - m_h[i] * xjui0 / 6.;
        ++j;
      }
    } // Could be optimized, e.g. with iterators, but executed only once
  }

  /**
   * @brief Interpolate and integrate given values.
   * @param y The knots y's
   * @param w The weights
   */
  template <typename T>
  T operator()(const T* y, const double* w);

private:
  struct Constants {
    double y0, y1, z0, z1;
  };

  /**
   * @brief The number of splines.
   */
  std::size_t m_size;

  /**
   * @brief The sample x's.
   */
  std::vector<double> m_u;

  /**
   * @brief The sample steps.
   */
  std::vector<double> m_h;

  /**
   * @brief The integration x's.
   */
  std::vector<double> m_x;

  /**
   * @brief The number of integration x's per spline.
   */
  std::vector<long> m_i;

  /**
   * @brief The pre-computed constants.
   */
  std::vector<Constants> m_k;
};

} // namespace Miller
} // namespace Phi

#endif // _PHIMILLER_LAMBDAINTEGRATOR_H
