// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef _PHIZERNIKE_CHONG_H
#define _PHIZERNIKE_CHONG_H

#include "LitlCore/Raster.h"
#include "PhiZernike/Indexing.h"

#include <cmath> // atan
#include <limits> // NaN
#include <vector>

namespace Phi {
namespace Zernike {

/**
 * @brief Radial polynomials at given radial distance.
 */
class ChongRadials {

private:
  /**
   * @brief Powers of r.
   */
  std::vector<double> powers;

  /**
   * @brief Polynomials with ANSI indexing.
   */
  std::vector<double> radials;

  /**
   * @brief Azimuthal degrees.
   */
  std::vector<long> ms;

  /**
   * @brief Current ANSI index.
   */
  std::size_t index;

  /**
   * @brief Max computed order.
   */
  long degree;

public:
  /**
   * @brief Constructor.
   */
  ChongRadials(double r2, long degrees) : powers {std::sqrt(r2), r2}, radials(), ms(), index(0), degree(-1) {
    upTo(degrees - 1);
  }

  const std::vector<double>& polynoms() const {
    return radials;
  }

  const std::vector<long>& azimuthals() const {
    return ms;
  }

  /**
   * @brief Get the Zernike radial polynomials up to order n in ANSI ordering.
   */
  const std::vector<double>& upTo(long n) {
    evalUpTo(n);
    return radials;
  }

private:
  /**
   * @brief Evaluate polynomials up to given order.
   */
  void evalUpTo(long n) {
    for (long p = degree + 1; p <= n; ++p) {
      evalOrder(p);
    }
    degree = n;
  }

  /**
   * @brief Evaluate polynomials at given order.
   * @details
   * Assuming lower orders have been evaluated already.
   */
  void evalOrder(long p) {
    long forward = radials.size();
    index = forward + p;
    radials.resize(index + 1);
    ms.resize(index + 1);
    for (long q = p; q >= 0; q -= 2, --index, ++forward) {
      ms[index] = q;
      const auto value = evalAt(p, q);
      radials[index] = value;
      if (q != 0) {
        ms[forward] = -q;
        radials[forward] = value;
      }
    }
  }

  /**
   * @brief Evaluate polynomial at given order and repetition.
   * @details
   * Assuming lower orders and higher repetitions have been evaluated already.
   */
  double evalAt(long p, long q) {
    if (p == 0) {
      return 1;
    }
    if (p == 1) {
      return pow(1);
    }
    if (p == 2 && q == 0) {
      return 2. * pow(2) - 1;
    }
    if (p == q) {
      return pow(p);
    }
    if (p - q == 2) {
      return p * at(p, p) - (p - 1) * at(p - 2, p - 2); // FIXME compute indices
    }
    const double h3 = -4 * (q - 2) * (q - 3) / ((p + q - 2) * (p - q + 4));
    const double h2 = h3 * (p + q) * (p - q + 2) / (4 * (q - 1)) + q - 2;
    const double h1 = q * (q - 1) - q * h2 + h3 * (p + q + 2) * (p - q) / 8;
    return h1 * radials[index + 2] + (h2 + h3 / pow(2)) * radials[index + 1];
  }

  /**
   * @brief p-th power of r.
   * @details Lazy evaluation.
   */
  double pow(long p) {
    while (powers.size() < std::size_t(p)) {
      powers.push_back(powers[powers.size() - 1] * powers[0]);
    }
    return powers[p - 1];
  }

  /**
   * @brief Value at given degrees.
   */
  double& at(long p, long q) {
    return radials[as<AnsiIndex>(ZernikeDegrees {p, q}).index];
  }
};

/**
 * @brief Zernike basis.
 */
Litl::Raster<double, 3> chongBasis(long side, long orders) {
  const ChongRadials zero(0, orders - 1);
  const auto& azimuths = zero.azimuthals();
  const long count = azimuths.size();
  Litl::Raster<double, 3> raster({count, side, side});
  const double center = .5 * side;
  const double normalization = 1. / (center * center);
  double* it = raster.data();
  for (long y = -center; y < center; ++y) {
    for (long x = -center; x < center; ++x) {
      const double r2 = normalization * (x * x + y * y);
      const double theta = std::atan2(y, x);
      if (r2 <= 1) {
        const ChongRadials radials(r2, orders);
        const auto& values = radials.polynoms();
        for (long i = 0; i < count; ++i, ++it) {
          const long m = azimuths[i];
          double sine = 1;
          if (m > 0) {
            sine = std::cos(theta * m);
          } else if (m < 0) {
            sine = std::sin(theta * m);
          }
          *it = sine * values[i];
        }
      } else {
        std::fill(it, it + count, 0);
        it += count;
      }
    }
  }
  return raster;
}

} // namespace Zernike
} // namespace Phi

#endif // _PHIZERNIKE_CHONG_H
