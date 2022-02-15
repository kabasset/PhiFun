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

#ifndef _PHIZERNIKE_INDEXING_H
#define _PHIZERNIKE_INDEXING_H

namespace Phi {
namespace Zernike {

/**
 * @brief Zernike original indexing (radial and azimuthal degrees).
 */
struct ZernikeDegrees {
  long radial, azimuthal;
};

/**
 * @brief OSA/ANSI index.
 */
struct AnsiIndex {
  long index;
};

/**
 * @brief Noll index.
 */
struct NollIndex {
  long index;
};

/**
 * @brief Convert indexing scheme.
 */
template <typename TTo, typename TFrom>
constexpr TTo as(const TFrom& index);

/**
 * @brief Zernike degrees to ANSI index.
 */
template <>
constexpr AnsiIndex as(const ZernikeDegrees& degrees) {
  return AnsiIndex {(degrees.radial * (degrees.radial + 2) + degrees.azimuthal) / 2};
}

/**
 * @brief Noll index to Zernike degrees.
 */
template <>
constexpr ZernikeDegrees as(const NollIndex& noll) {
  long radial = 0;
  long i = noll.index - 1;
  while (i > radial) {
    ++radial;
    i -= radial;
  }
  long azimuthal = (radial % 2) + 2 * long((i + ((radial + 1) % 2)) / 2);
  if (noll.index % 2 == 0) {
    return ZernikeDegrees {radial, azimuthal};
  } else {
    return ZernikeDegrees {radial, -azimuthal};
  }
}

/**
 * @brief Noll index to ANSI index.
 */
template <>
constexpr AnsiIndex as(const NollIndex& noll) {
  return as<AnsiIndex>(as<ZernikeDegrees>(noll));
}

} // namespace Zernike
} // namespace Phi

#endif // _PHIZERNIKE_INDEXING_H
