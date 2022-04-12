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

template <typename T>
struct Kernel1D {

  Kernel1D(const std::vector<T>& values, long anchor) :
      backward(anchor), forward(values.size() - 1 - backward), center(&values[backward]), bias() {}

  long backward;
  long forward;
  const T* center;
  T bias;

  long size() const {
    return backward + forward + 1;
  }

  const T* begin() const {
    return center - backward;
  }

  const T* end() const {
    return center + forward + 1;
  }
};

template <typename T>
struct Sampling1D : public std::iterator<std::input_iterator_tag, T> {

  Sampling1D(long size, T* data) :
      m_size(size), m_from(0), m_to(m_size - 1), m_step(1), m_stride(1), m_data(data), m_it(m_data + m_from) {}

  std::size_t size() const {
    return m_size;
  }

  long count() const {
    return (m_to - m_from + 1) / (m_step * m_stride) * (m_step * m_stride);
  }

  long from() const {
    return m_from;
  }

  void from(long value) {
    m_from = value;
  }

  long to() const {
    return m_to;
  }

  void to(long value) {
    m_to = value;
  }

  long step() const {
    return m_step;
  }

  void step(long value) {
    m_step = value;
  }

  long stride() const {
    return m_stride;
  }

  void stride(long value) {
    m_stride = value;
  }

  Sampling1D<const T> begin() const {
    Sampling1D<const T> b(m_size, m_data);
    b.m_from = m_from;
    b.m_to = m_to;
    b.m_step = m_step;
    b.m_stride = m_stride;
    return b;
  }

  Sampling1D<const T> end() const {
    auto e = begin();
    e.m_it = (m_data + m_from) + count();
    return e;
  }
  Sampling1D<T> begin() {
    auto b = *this;
    b.m_it = m_data + m_from;
    return b;
  }

  Sampling1D<T> end() {
    auto e = *this;
    e.m_it = (m_data + m_from) + count();
    return e;
  }

  T& operator*() {
    return *m_it;
  }

  T* operator->() {
    return m_it;
  }

  Sampling1D<T>& operator++() {
    m_it += m_step * m_stride;
    return *this;
  }

  Sampling1D<T> operator++(int) {
    auto res = *this;
    ++res;
  }

  Sampling1D<T>& operator+=(long n) {
    m_it += m_step * m_stride * n;
    return *this;
  }

  bool operator==(const Sampling1D<T>& rhs) const {
    return m_it == rhs.m_it;
  }

  bool operator!=(const Sampling1D& rhs) const {
    return not(*this == rhs);
  }

private:
  long m_size;
  long m_from;
  long m_to;
  long m_step;
  long m_stride;
  T* m_data;
  T* m_it;
};

template <typename TSampling, typename TKernel, typename USampling>
void convolve1DZero(const TSampling& in, const TKernel& kernel, USampling& out) {
  auto inIt = in;
  inIt.step(1);
  auto outIt = out;
  // outIt.step(1);
  long i = in.from();
  for (; i < kernel.backward; i += in.step(), outIt += out.step()) {
    *outIt = std::inner_product(kernel.center - i, kernel.end(), inIt, kernel.bias);
  }
  for (; i < in.size() - kernel.forward; i += in.step(), inIt += in.step(), outIt += out.step()) {
    *outIt = std::inner_product(kernel.begin(), kernel.end(), inIt, kernel.bias);
  }
  for (; i <= in.to(); i += in.step(), inIt += in.step(), outIt += out.step()) {
    *outIt = std::inner_product(kernel.begin(), kernel.center + (in.size() - 1 - i), inIt, kernel.bias);
  }
}

template <typename TSampling, typename TKernel, typename USampling>
void convolveXyZero(const TSampling& in, const TKernel& kernel, USampling& out) {}

} // namespace Image2D
} // namespace Phi

#endif // _PHIBOX_IMAGEPROCESSING_H
