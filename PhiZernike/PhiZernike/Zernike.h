// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef _PHIZERNIKE_ZERNIKE_H
#define _PHIZERNIKE_ZERNIKE_H

#include "EleFitsData/Raster.h"
#include "PhiZernike/Indexing.h"

#include <algorithm> // copy_n
#include <array>
#include <limits> // NaN
#include <utility> // index_sequence

namespace Phi {
namespace Zernike {

/**
 * @brief Zernike polynomials calculator at given (x, y).
 * @warning
 * The computation is limited to `JCount` values.
 */
class XySeries {
public:
  static constexpr long JCount = 45;

public:
  /**
   * @brief Constructor.
   * @details
   * Computes intermediate variables only.
   */
  XySeries(double u, double v, double blank = std::numeric_limits<double>::quiet_NaN()) :
      x1(v), // inverted
      x2(x1 * x1), x3(x1 * x2), x4(x1 * x3), x5(x1 * x4), x6(x1 * x5), x7(x1 * x6), x8(x1 * x7), x9(x1 * x8),
      y1(u), // inverted
      y2(y1 * y1), y3(y1 * y2), y4(y1 * y3), y5(y1 * y4), y6(y1 * y5), y7(y1 * y6), y8(y1 * y7), y9(y1 * y8),
      nan(blank) {}

  /**
   * @brief Evaluate at given indices.
   * @details
   * \code
   * series.at<NollIndex, 4, 11>();
   * \endcode
   */
  template <typename TIndex, long... Js>
  std::array<double, sizeof...(Js)> at(std::integer_sequence<long, Js...> = {}) {
    return {eval<TIndex, Js>()...};
  }

  /**
   * @brief Evaluate up to given index.
   */
  template <typename TIndex, long JMax>
  std::array<double, JMax + 1> upTo() {
    return at<TIndex>(std::make_integer_sequence<long, JMax + 1>());
  }

  /**
   * @brief Evaluate between two included indices.
   */
  template <typename TIndex, long From, long To>
  std::array<double, To - From + 1> fromTo() {
    return at<TIndex>(offset<From>(std::make_integer_sequence<long, To - From + 1>()));
  }

  /**
   * @brief Evaluate between two included indices.
   */
  template <typename TIndex>
  std::vector<double> fromTo(long from, long to) {
    std::vector<double> v(to - from + 1);
    fromTo<TIndex>(v);
    return v;
  }

  /**
   * @brief Evaluate between two included indices into some data array.
   */
  template <typename TIndex>
  void fromTo(long from, long to, double* data) {
    const auto a = upTo<TIndex, JCount - 1>();
    const auto n = to - from + 1;
    std::copy_n(a.begin() + from, n, data);
  }

private:
  /**
   * @brief Evaluate at given index.
   * @details
   * \code
   * series.eval<NollIndex, 11>();
   * \endcode
   */
  template <typename TIndex, long J>
  double eval() const {
    return ansi<as<AnsiIndex>(TIndex {J}).index>();
  }

  /**
   * @brief Evaluate at given ANSI index.
   * @details
   * \code
   * series.eval<11>();
   * \endcode
   */
  template <long J>
  double ansi() const;

  /**
   * @brief Offset a `std::integer sequence`.
   */
  template <long Offset, long... Js>
  std::index_sequence<(Offset + Js)...> offset(std::index_sequence<Js...>) {
    return {};
  }

  /**
   * @brief Powers of x.
   */
  double x1, x2, x3, x4, x5, x6, x7, x8, x9;

  /**
   * @brief Powers of y.
   */
  double y1, y2, y3, y4, y5, y6, y7, y8, y9;

  /**
   * @brief Value for |(x, y)| > 1.
   */
  double nan;

  std::array<double, JCount> series;
};

/**
 * @brief Zernike basis.
 */
Euclid::Fits::VecRaster<double, 3> basis(long side, long count = XySeries::JCount) {
  Euclid::Fits::VecRaster<double, 3> raster({count, side, side});
  const double center = .5 * side;
  const double normalization = 1. / center;
  double* it = raster.data();
  for (long y = -center; y < center; ++y) {
    for (long x = -center; x < center; ++x, it += count) {
      const double u = normalization * x;
      const double v = normalization * y;
      if (u * u + v * v > 1) {
        std::fill(it, it + count, 0);
      } else {
        XySeries series(u, v);
        series.fromTo<NollIndex>(0, count - 1, it); // FIXME allow any ordering
      }
    }
  }
  return raster;
}

} // namespace Zernike
} // namespace Phi

#include "PhiZernike/impl/Zernike.hpp"

#endif // _PHIZERNIKE_ZERNIKE_H
