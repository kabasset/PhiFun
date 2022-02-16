/**
 * @copyright (C) 2012-2020 Euclid Science Ground Segment
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 3.0 of the License, or (at your option)
 * any later version.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 *
 */

#ifndef _PHIZERNIKE_IMPL_ZERNIKE_HPP
#define _PHIZERNIKE_IMPL_ZERNIKE_HPP

namespace Phi {
namespace Zernike {

// Precompute first polynomials

#define DEF_ZERNIKE(J, expr) \
  template <> \
  double XySeries::ansi<J>() const { \
    if (x2 + y2 > 1) { \
      return nan; \
    } \
    return expr; \
  }

// https://jumar.lowell.edu/confluence/download/attachments/20546539/2011.JMOp.58.545L.pdf?version=2&modificationDate=1610084878000&api=v2
// See (Noll ordering): https://arxiv.org/pdf/0809.2368.pdf

// 0
DEF_ZERNIKE(0, 1)
// 1
DEF_ZERNIKE(1, x1)
DEF_ZERNIKE(2, y1)
// 2
DEF_ZERNIKE(3, 2 * x1 * y1)
DEF_ZERNIKE(4, -1 + 2 * x2 + 2 * y2)
DEF_ZERNIKE(5, -x2 + y2)
// 3
DEF_ZERNIKE(6, -x3 + 3 * x1 * y2)
DEF_ZERNIKE(7, -2 * x1 + 3 * x3 + 3 * x1 * y2)
DEF_ZERNIKE(8, -2 * y1 + 3 * y3 + 3 * x2 * y1)
DEF_ZERNIKE(9, y3 - 3 * x2 * y1)
// 4
DEF_ZERNIKE(10, -4 * x3 * y1 + 4 * x1 * y3)
DEF_ZERNIKE(11, -6 * x1 * y1 + 8 * x3 * y1 + 8 * x1 * y3)
DEF_ZERNIKE(12, 1 - 6 * x2 - 6 * y2 + 6 * x4 + 12 * x2 * y2 + 6 * y4)
DEF_ZERNIKE(13, 3 * x2 - 3 * y2 - 4 * x4 + 4 * y4)
DEF_ZERNIKE(14, x4 - 6 * x2 * y2 + y4)
// 5
DEF_ZERNIKE(15, x5 - 10 * x3 * y2 + 5 * x1 * y4)
DEF_ZERNIKE(16, 4 * x3 - 12 * x1 * y2 - 5 * x5 + 10 * x3 * y2 + 15 * x1 * y4)
DEF_ZERNIKE(17, 3 * x1 - 12 * x3 - 12 * x1 * y2 + 10 * x5 + 20 * x3 * y2 + 10 * x1 * y4)
DEF_ZERNIKE(18, 3 * y1 - 12 * y3 - 12 * x2 * y1 + 10 * y5 + 20 * x2 * y3 - 15 * x4 * y1)
DEF_ZERNIKE(19, -4 * y3 + 12 * x2 * y1 + 5 * y5 - 10 * x2 * y3 - 15 * x4 * y1)
DEF_ZERNIKE(20, y5 - 10 * x2 * y3 + 5 * x4 * y1)
// 6
DEF_ZERNIKE(21, 6 * x5 * y1 - 20 * x3 * y3 + 6 * x1 * y5)
DEF_ZERNIKE(22, 20 * x3 * y1 - 20 * x1 * y3 - 24 * x5 * y1 + 24 * x1 * y5)
DEF_ZERNIKE(23, 12 * x1 * y1 - 40 * x3 * y1 - 40 * x1 * y3 + 30 * x5 * y1 + 60 * x3 * y3 - 30 * x1 * y5)
DEF_ZERNIKE(
    24,
    -1 + 12 * x2 + 12 * y2 - 30 * x4 - 60 * x2 * y2 - 30 * y4 + 20 * x6 + 60 * x4 * y2 + 60 * x2 * y4 + 20 * y6)
DEF_ZERNIKE(25, -6 * x2 + 6 * y2 + 20 * x4 - 20 * y4 - 15 * x6 - 15 * x4 * y2 + 15 * x2 * y4 + 15 * y6)
DEF_ZERNIKE(26, -5 * x4 + 30 * x2 * y2 - 5 * y4 + 6 * x6 - 30 * x4 * y2 - 30 * x2 * y4 + 6 * y6)
DEF_ZERNIKE(27, -x6 + 15 * x4 * y2 - 15 * x2 * y4 + y6)
// 7
DEF_ZERNIKE(28, -x7 + 21 * x5 * y2 - 35 * x3 * y4 + 7 * x1 * y6)
DEF_ZERNIKE(29, -6 * x5 + 60 * x3 * y2 - 30 * x1 * y4 + 7 * x7 - 63 * x5 * y2 - 35 * x3 * y4 + 35 * x1 * y6)
DEF_ZERNIKE(
    30,
    -10 * x3 + 30 * x1 * y2 + 30 * x5 - 60 * x3 * y2 - 90 * x1 * y4 - 21 * x7 + 21 * x5 * y2 + 105 * x3 * y4 +
        63 * x1 * y6)
DEF_ZERNIKE(
    31,
    -4 * x1 + 30 * x3 + 30 * x1 * y2 - 60 * x5 - 120 * x3 * y2 - 60 * x1 * y4 + 35 * x7 + 105 * x5 * y2 +
        105 * x3 * y4 + 35 * x1 * y6)
DEF_ZERNIKE(
    32,
    -4 * y1 + 30 * y3 + 30 * x2 * y1 - 60 * y5 - 120 * x2 * y3 - 60 * x4 * y1 + 35 * y7 + 105 * x2 * y5 +
        105 * x4 * y3 + 35 * x6 * y1)
DEF_ZERNIKE(
    33,
    30 * y3 - 30 * x2 * y1 - 30 * y5 + 60 * x2 * y3 + 90 * x4 * y1 + 21 * y7 - 21 * x2 * y5 - 105 * x4 * y3 +
        63 * x6 * y1)
DEF_ZERNIKE(34, -6 * y5 + 60 * x2 * y3 - 30 * x4 * y1 + 7 * y7 - 63 * x2 * y5 - 35 * x4 * y3 + 35 * x6 * y1)
DEF_ZERNIKE(35, x7 - 21 * x2 * y5 + 35 * x4 * y3 - 7 * x6 * y1)
// 8
DEF_ZERNIKE(37, -42 * x5 * y1 + 140 * x3 * y3 - 42 * x1 * y5 + 48 * x7 * y1 - 112 * x5 * y3 + 48 * x1 * y7)
DEF_ZERNIKE(
    38,
    -60 * x3 * y1 + 60 * x1 * y3 + 168 * x5 * y1 - 168 * x1 * y5 - 112 * x7 * y1 - 112 * x5 * y3 + 112 * x3 * y5 +
        112 * x1 * y7)
DEF_ZERNIKE(
    39,
    -20 * x1 * y1 + 120 * x3 * y1 + 120 * x1 * y3 - 210 * x5 * y1 - 420 * x3 * y3 - 210 * x1 * y5 - 112 * x7 * y1 +
        336 * x5 * y3 + 336 * x3 * y5 + 112 * x1 * y7)
DEF_ZERNIKE(
    40,
    1 - 20 * x2 - 20 * y2 + 90 * x4 + 180 * x2 * y2 + 90 * y4 - 140 * x6 - 420 * x4 * y2 - 420 * x2 * y4 - 140 * y6 +
        70 * x8 + 280 * x6 * y2 + 420 * x4 * y4 + 280 * x2 * y6 + 70 * y8)
DEF_ZERNIKE(
    41,
    10 * x2 - 10 * y2 - 60 * x4 + 105 * x4 * y2 - 105 * x2 * y4 + 60 * y4 + 105 * x6 - 105 * y6 - 56 * x8 -
        112 * x6 * y2 + 112 * x2 * y6 + 56 * y8)
DEF_ZERNIKE(
    42,
    15 * x4 - 90 * x2 * y2 + 15 * y4 - 42 * x6 + 210 * x4 * y2 + 210 * x2 * y4 - 42 * y6 + 28 * x8 - 112 * x6 * y2 -
        280 * x4 * y4 - 112 * x2 * y6 + 28 * y8)
DEF_ZERNIKE(43, 7 * x6 - 105 * x4 * y2 + 105 * x2 * y4 - 7 * y6 - 8 * x8 + 112 * x6 * y2 - 112 * x2 * y6 + 8 * y8)
DEF_ZERNIKE(44, x8 - 28 * x6 * y2 + 70 * x4 * y4 - 28 * x2 * y6 + y8)

#undef DEF_ZERNIKE

} // namespace Zernike
} // namespace Phi

#endif // _PHIZERNIKE_IMPL_ZERNIKE_HPP
