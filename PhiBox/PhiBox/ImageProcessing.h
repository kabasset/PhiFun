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
U bilinearNearest(const TRaster& raster, double x, double y, TLongs... coords) {
  const auto width = raster.shape()[0];
  const auto height = raster.shape()[1];
  const auto left = clamp<long>(x, 0, width - 1);
  const auto bottom = clamp<long>(y, 0, height - 1);
  const auto* bl = &raster[{left, bottom, coords...}];
  if (x == left && y == bottom) { // No interpolation
    return *bl;
  }
  const double dx = x - left;
  const double dy = y - bottom;
  const auto* br = x >= width - 1 ? bl : bl + 1;
  const auto* tl = y >= height - 1 ? bl : bl + width;
  const auto* tr = x >= width - 1 ? tl : tl + 1;
  return *bl + (*br - *bl + (*bl + *tr - *br - *tl) * dy) * dx + (*tl - *bl) * dy;
}

template <typename U, typename TRaster, typename... TLongs>
U bilinearZero(const TRaster& raster, double x, double y, TLongs... coords) {
  const auto width = raster.shape()[0];
  const auto height = raster.shape()[1];
  if (x < 0 || x > width - 1 || y < 0 || y > height - 1) {
    return 0;
  }
  const auto left = long(x);
  const auto bottom = long(y);
  const auto* bl = &raster[{left, bottom, coords...}];
  if (x == left && y == bottom) { // No interpolation
    return *bl;
  }
  const double dx = x - left;
  const double dy = y - bottom;
  const auto* br = bl + 1;
  const auto* tl = bl + width;
  const auto* tr = tl + 1;
  return *bl + (*br - *bl + (*bl + *tr - *br - *tl) * dy) * dx + (*tl - *bl) * dy;
}

} // namespace Image2D
} // namespace Phi

#endif // _PHIBOX_IMAGEPROCESSING_H
