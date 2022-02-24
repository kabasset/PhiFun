// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "EleFitsData/Raster.h"

#ifndef _PHIBOX_IMAGEPROCESSING_H
#define _PHIBOX_IMAGEPROCESSING_H

namespace Phi {
namespace Image2D {

template <typename T>
T clamp(T in, T min, T max) {
  return std::max(min, std::min(in, max));
}

template <typename U, typename TRaster, typename... TLongs>
U bilinear(const TRaster& raster, double x, double y, TLongs... coords) {
  using T = typename TRaster::Value;
  const auto width = raster.shape()[0];
  const auto height = raster.shape()[1];
  const auto left = clamp<long>(x, 0, width - 1);
  const auto bottom = clamp<long>(y, 0, height - 1);
  const T* bl = &raster[{left, bottom, coords...}];
  if (x == left && y == bottom) { // No interpolation
    return *bl;
  }
  const double dx = x - left;
  const double dy = y - bottom;
  const T* br = x >= width - 1 ? bl : bl + 1;
  const T* tl = y >= height - 1 ? bl : bl + width;
  const T* tr = x >= width - 1 ? tl : tl + 1;
  return *bl + (*br - *bl + (*bl + *tr - *br - *tl) * dy) * dx + (*tl - *bl) * dy;
}

} // namespace Image2D
} // namespace Phi

#endif // _PHIBOX_IMAGEPROCESSING_H
