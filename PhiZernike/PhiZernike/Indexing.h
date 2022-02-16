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
 * @brief Fringe/University of Arizona index.
 */
struct FringeIndex {
  long index;
};

/**
 * @brief Wyant index.
 */
struct WyantIndex {
  long index;
};

/**
 * @brief Convert indexing scheme from or to Zernike degrees.
 * @details
 * \code
 * const auto noll = as<NollIndex>(ZernikeDegrees {3, -1});
 * \endcode
 */
template <typename TTo, typename TFrom>
constexpr TTo as(const TFrom& index);

} // namespace Zernike
} // namespace Phi

#include "PhiZernike/impl/Indexing.hpp"

#endif // _PHIZERNIKE_INDEXING_H
