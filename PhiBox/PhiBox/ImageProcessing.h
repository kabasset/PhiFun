// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "EleFitsData/Raster.h"

#ifndef _PHIBOX_IMAGEPROCESSING_H
#define _PHIBOX_IMAGEPROCESSING_H

namespace Phi {
namespace Box {

template <typename U, typename TRaster>
U bilinear(const TRaster& raster, double x, double y) {
  using T = typename TRaster::Value;
  const long left = x;
  const long bottom = y;
  const T* bl = &raster[{left, bottom}];
  if (x == left && y == bottom) { // No interpolation
    return *bl;
  }
  const double dx = x - left;
  const double dy = y - bottom;
  const T* br = bl + 1;
  const T* tl = bl + raster.template length<0>();
  const T* tr = tl + 1;
  return *bl + (*br - *bl + (*bl + *tr - *br - *tl) * dy) * dx + (*tl - *bl) * dy;
}

} // namespace Box
} // namespace Phi

#endif // _PHIBOX_IMAGEPROCESSING_H
