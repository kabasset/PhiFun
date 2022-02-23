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
   * @brief Compute the second derivatives at knots.
   */
  template <typename T>
  std::vector<T> knotZ(const T* y) {
    const auto m_h0 = m_h[0]; // FIXME as members
    const auto m_h1 = m_h[1];
    const auto m_h02 = m_h0 * m_h0;
    const auto m_h12 = m_h1 * m_h1;
    const auto m_h01 = m_h0 * m_h1;
    const auto m_hm0 = m_h[m_size - 1];
    const auto m_hm1 = m_h[m_size - 2];
    const auto m_hm02 = m_hm0 * m_hm0;
    const auto m_hm12 = m_hm1 * m_hm1;
    const auto m_hm01 = m_hm0 * m_hm1;

    std::vector<T> s(m_size);
    std::vector<T> z(m_size);
    s[0] = (y[1] - y[0]) / m_h0; // Because next loop starts at 1
    for (std::size_t i = 1; i < m_size - 1; ++i) {
      s[i] = (y[i + 1] - y[i]) / m_h[i];
      z[i] = (y[i + 1] - y[i]) / m_h[i] - (y[i] - y[i - 1]) / m_h[i - 1];
    } // FIXME optimize with iterators

    // Not-a-knot at 0
    const auto s0 = s[0];
    const auto s1 = s[1];
    const auto z1 = z[1];
    z[0] = (6. * m_h1 * (s1 - s0) - (2 * m_h02 + 3 * m_h01 + m_h12) * z1) / (m_h02 - m_h12);

    // Not-a-knot at m_size - 1
    const auto sm0 = s[m_size - 1];
    const auto sm1 = s[m_size - 2];
    const auto zm1 = z[m_size - 2];
    z[m_size - 1] = (6. * m_hm0 * (sm0 - sm1) - (m_hm02 - m_hm12) * zm1) / (2. * m_hm02 + 3 * m_hm01 + m_hm12);

    return z;
  }

  /**
   * @brief Interpolate.
   * @param y The values at knots
   * @param z The second derivatives at knots
   */
  template <typename T>
  std::vector<T> interpolate(const T* y, const T* z) {
    std::size_t j = 0;
    for (std::size_t i = 0; i < m_size; ++i) {
      for (std::size_t n = 0; n < m_i[i]; ++n) {
        // FIXME constant-based formulae of Edo
        ++j;
      }
    }
  }

  /**
   * @brief Interpolate and integrate given values.
   * @param y The values at knots
   * @param z The second derivatives at knots
   * @param w The weights
   */
  template <typename T>
  T integrate(const T* y, const T* z, const double* w);

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
