// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

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
