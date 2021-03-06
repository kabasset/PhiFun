// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef _PHIBOX_SPLINEINTEGRATOR_H
#define _PHIBOX_SPLINEINTEGRATOR_H

#include <complex>
#include <vector>

namespace Phi {
namespace Spline {

/**
 * @brief Generate evenly spaced samples over a closed interval.
 */
std::vector<double> linspace(double min, double max, std::size_t count) {
  std::vector<double> xs(count);
  const double step = (max - min) / (count - 1);
  double value = min;
  for (auto& x : xs) {
    x = value;
    value += step;
  }
  return xs;
}

/**
 * @brief Spline interpolation and integration, inspired by a proposal by Lance Miller.
 * @details
 * Although this class is generic, it was initially developped for TF integration.
 * 
 * The monochromatic TFs are interpolated along lambdas at a few nanometers,
 * and then integrated along lambdas.
 * Lance Miller proposes to perform both operations at once,
 * i.e. for each position, compute the weight of interpolating splines and perform the summation.
 * The implementation is very different from that of Lance Miller,
 * but the intent is similar and the results are the same.
 * 
 * For example, go directly from 40 monochromatic TFs to 1 broadband TF
 * instead of computing 200 interpolated TFs in-between.
 */
class SplineIntegrator {

public:
  /**
   * @brief Constructor.
   * @param u The input knots x's (strictly increasing)
   * @param x The integration x's (strictly increasing and in u's interval)
   * @details
   * Precompute spline coefficients.
   */
  SplineIntegrator(std::vector<double> u, std::vector<double> x) :
      m_size(u.size() - 1), m_u(std::move(u)), m_h(m_size), m_g(m_size), m_x(std::move(x)), m_n(m_size),
      m_k(m_x.size()) {
    // FIXME assert m_size > 3
    std::size_t j = 0;
    for (std::size_t i = 0; i < m_size; ++i) {
      const double hi = m_u[i + 1] - m_u[i];
      // FIXME assert hi > 0
      m_h[i] = hi;
      m_g[i] = 1. / hi;
      while (j < m_k.size() && m_x[j] < m_u[i + 1]) {
        const double left = m_x[j] - m_u[i];
        const double right = m_u[i + 1] - m_x[j];
        m_k[j].z1 = left * left * left / (6. * hi) - hi * left / 6.;
        m_k[j].z0 = right * right * right / (6. * hi) - hi * right / 6.;
        m_k[j].y1 = left / hi;
        m_k[j].y0 = right / hi;
        ++m_n[i];
        ++j;
      }
    } // Could be optimized, e.g. with iterators, but executed only once and more readable this way

    m_h0 = m_h[0];
    m_h1 = m_h[1];
    m_h02 = m_h0 * m_h0;
    m_h12 = m_h1 * m_h1;
    m_h01 = m_h0 * m_h1;
    m_hm0 = m_h[m_size - 1];
    m_hm1 = m_h[m_size - 2];
    m_hm02 = m_hm0 * m_hm0;
    m_hm12 = m_hm1 * m_hm1;
    m_hm01 = m_hm0 * m_hm1;
    m_g0 = m_g[0];
    m_zm0Factor = 1. / (2. * m_hm0 * m_hm0 + 3 * m_hm0 * m_hm1 + m_hm1 * m_hm1);
  }

  /**
   * @brief Get the knot x's.
   */
  const std::vector<double>& knotX() const {
    return m_u;
  }

  /**
   * @brief Get the interpolation x's.
   */
  const std::vector<double>& interpolationX() const {
    return m_x;
  }

  /**
   * @brief Compute the second derivatives at knots.
   */
  template <typename T>
  std::vector<T> knotZ(const T* y) const {
    std::vector<T> s(m_size);
    std::vector<T> z(m_size + 1);
    s[0] = (y[1] - y[0]) * m_g0; // Because next loop starts at 1
    auto* yIt = y + 1;
    auto* gIt = m_g.data() + 1;
    auto* sIt = s.data() + 1;
    auto* zIt = z.data() + 1;
    for (auto i = m_size - 1; i--; ++yIt, ++gIt, ++sIt, ++zIt) {
      // s[i] = (y[i + 1] - y[i]) / m_h[i];
      *sIt = (*(yIt + 1) - *yIt) * *gIt;
      // z[i] = (y[i + 1] - y[i]) / m_h[i] - (y[i] - y[i - 1]) / m_h[i - 1];
      *zIt = *sIt - *(sIt - 1);
    }

    // Not-a-knot at 0
    const auto s0 = s[0];
    const auto s1 = s[1];
    const auto z1 = z[1];
    const auto z2 = z[2];
    z[0] = (6. * (s1 - s0) - 2. * (m_h0 + m_h1) * z1 - m_h1 * z2) * m_g0;

    // Not-a-knot at m_size - 1
    // FIXME Should it be at m_size?
    const auto sm0 = s[m_size - 1];
    const auto sm1 = s[m_size - 2];
    const auto zm1 = z[m_size - 1];
    z[m_size - 1] = (6. * m_hm0 * (sm0 - sm1) - (m_hm02 - m_hm12) * zm1) * m_zm0Factor;

    return z;
  }

  /**
   * @brief Interpolate.
   * @param y The values at knots
   * @param z The second derivatives at knots
   */
  template <typename T>
  std::vector<T> interpolate(const T* y, const T* z) const {
    std::vector<T> s(m_x.size());
    std::size_t j = 0;
    for (std::size_t i = 0; i < m_size; ++i) {
      for (long n = 0; n < m_n[i]; ++n) {
        s[j] = z[i + 1] * m_k[j].z1 + z[i] * m_k[j].z0 + y[i + 1] * m_k[j].y1 + y[i] * m_k[j].y0;
        ++j;
      }
    }
    return s;
  }

  /**
   * @brief Interpolate and integrate.
   * @param y The values at knots
   * @param z The second derivatives at knots
   * @param w The weights
   */
  template <typename T>
  T integrate(const T* y, const T* z, const double* w) const {
    T sum {};
    auto* yIt = y;
    auto* zIt = z;
    auto* wIt = w;
    auto kIt = m_k.begin();
    for (auto i = m_size; i--; ++yIt, ++zIt) {
      for (auto n = m_n[i]; n--; ++wIt, ++kIt) {
        // sum += w[j] * (z[i + 1] * m_k[i].z1 + z[i] * m_k[i].z0 + y[i + 1] * m_k[i].y1 + y[i] * m_k[i].y0);
        sum += *wIt * (*(zIt + 1) * kIt->z1 + *zIt * kIt->z0 + *(yIt + 1) * kIt->y1 + *yIt * kIt->y0);
      }
    }
    return sum;
  }

private:
  /**
   * @brief The coefficients of a polynom.
   */
  struct Coefficients {
    double z1; ///< Coefficient of z[i + 1]
    double z0; ///< Coefficient of z[i]
    double y1; ///< Coefficient of y[i + 1]
    double y0; ///< Coefficient of y[i]
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
  double m_h0; ///< m_h[0]
  double m_h1; ///< m_h[1]
  double m_h02; ///< m_h0 * m_h0
  double m_h12; ///< m_h1 * m_h1
  double m_h01; ///< m_h0 * m_h1
  double m_hm0; ///< m_h[m_size - 1]
  double m_hm1; ///< m_h[m_size - 2]
  double m_hm02; ///< m_hm0 * m_hm0
  double m_hm12; ///< m_hm1 * m_hm1
  double m_hm01; ///< m_hm0 * m_hm1
  std::vector<double> m_g; ///< Inverses of m_h[i]'s
  double m_g0; ///< m_g[0]
  double m_zm0Factor; ///< Normalization factor for z[m_size - 1]

  /**
   * @brief The integration x's.
   */
  std::vector<double> m_x;

  /**
   * @brief The number of integration x's per spline.
   */
  std::vector<long> m_n;

  /**
   * @brief The pre-computed constants.
   */
  std::vector<Coefficients> m_k;
};

} // namespace Spline
} // namespace Phi

#endif // _PHIBOX_SPLINEINTEGRATOR_H
