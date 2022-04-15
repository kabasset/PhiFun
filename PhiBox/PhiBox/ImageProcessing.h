// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "EleFitsData/Raster.h"

#ifndef _PHIBOX_IMAGEPROCESSING_H
#define _PHIBOX_IMAGEPROCESSING_H

namespace Phi {
namespace Image2D {

using Euclid::Fits::Position;

using Euclid::Fits::Region;

using Euclid::Fits::VecRaster;

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
struct Sampling1D {

  template <typename TSampling>
  class Iterator : public std::iterator<std::input_iterator_tag, T> {

  public:
    using Value = typename TSampling::Value;

    Iterator(TSampling& sampling, long index) :
        m_sampling(sampling), m_it(sampling.data() + index * sampling.stride()) {}

    Value& operator*() {
      return *m_it;
    }

    Value* operator->() {
      return m_it;
    }

    Iterator& operator++() {
      return *this += 1;
    }

    Iterator operator++(int) {
      auto res = *this;
      ++res;
    }

    Iterator& operator+=(long n) {
      m_it += m_sampling.step() * m_sampling.stride() * n;
      return *this;
    }

    Iterator& operator-=(long n) {
      *this += -n;
      return *this;
    }

    bool operator==(const Iterator& rhs) const {
      return m_it == rhs.m_it;
    }

    bool operator!=(const Iterator& rhs) const {
      return not(*this == rhs);
    }

    Iterator& operator=(Value* it) {
      m_it = it;
      return *this;
    }

  private:
    TSampling& m_sampling;
    Value* m_it;
  };

  using Value = T;

  Sampling1D(long size, T* data) : m_size(size), m_from(0), m_to(m_size - 1), m_step(1), m_stride(1), m_data(data) {}

  std::size_t size() const {
    return m_size;
  }

  const T* data() const {
    return m_data;
  }

  T* data() {
    return m_data;
  }

  Sampling1D& data(T* d) {
    m_data = d;
    return *this;
  }

  long count() const {
    return (m_to - m_from + 1) / m_step;
  }

  long from() const {
    return m_from;
  }

  Sampling1D& from(long value) {
    m_from = value;
    return *this;
  }

  long to() const {
    return m_to;
  }

  Sampling1D& to(long value) {
    m_to = value;
    return *this;
  }

  long step() const {
    return m_step;
  }

  Sampling1D& step(long value) {
    m_step = value;
    return *this;
  }

  long stride() const {
    return m_stride;
  }

  Sampling1D& stride(long value) {
    m_stride = value;
    return *this;
  }

  Iterator<const Sampling1D<const T>> begin() const {
    return {*this, m_from};
  }

  Iterator<const Sampling1D<const T>> end() const {
    return {*this, m_to + 1};
  }

  Iterator<Sampling1D<T>> begin() {
    return {*this, m_from};
  }

  Iterator<Sampling1D<T>> end() {
    return {*this, m_to + 1};
  }

private:
  long m_size;
  long m_from;
  long m_to;
  long m_step;
  long m_stride;
  T* m_data;
};

template <typename TSampling, typename TKernel, typename USampling>
void convolve1DZero(const TSampling& in, const TKernel& kernel, USampling& out) {
  auto inUnit = in;
  inUnit.step(1); // For inner_product
  auto inIt = inUnit.begin();
  auto inMinIt = inIt;
  inMinIt -= in.from();
  inIt -= kernel.backward;
  auto outIt = out.begin();
  long i = in.from();
  for (; i < kernel.backward; i += in.step(), inIt += in.step(), ++outIt) {
    *outIt = std::inner_product(kernel.center - i, kernel.end(), inMinIt, kernel.bias);
  }
  for (; i < in.size() - kernel.forward; i += in.step(), inIt += in.step(), ++outIt) {
    *outIt = std::inner_product(kernel.begin(), kernel.end(), inIt, kernel.bias);
  }
  for (; i <= in.to(); i += in.step(), inIt += in.step(), ++outIt) {
    *outIt = std::inner_product(kernel.begin(), kernel.center + (in.size() - i), inIt, kernel.bias);
  }
}

template <typename TRasterIn, typename TKernel, typename TRasterOut>
void convolveXyZero(
    const TRasterIn& in,
    const Region<2>& region,
    const Position<2>& step,
    const TKernel& kernel,
    TRasterOut& out) { // FIXME support nd-rasters for batch processing

  // Set sampling
  const auto xFrom = region.front[0];
  const auto yFrom = std::max(0L, region.front[1] - kernel.backward);
  const auto xTo = region.back[0];
  const auto yTo = std::min(in.template length<1>() - 1, region.back[1] + kernel.forward);
  const auto xStep = step[0];
  const auto yStep = step[1];
  printf("%i-%i:%i - %i-%i:%i\n", xFrom, xTo, xStep, yFrom, yTo, yStep);

  // Convolve along x-axis
  const long xConvolvedWidth = (xTo - xFrom + 1) / xStep;
  const long xConvolvedHeight = (yTo - yFrom + 1); // / yStep;
  VecRaster<typename TRasterOut::Value, TRasterOut::Dim> xConvolved({xConvolvedWidth, xConvolvedHeight});
  Sampling1D<const typename TRasterIn::Value> inSampling(region.shape()[0], &in[{0, yFrom}]);
  inSampling.from(xFrom).to(xTo).step(xStep);
  printf("inSampling: %i-%i:%i:%i\n", inSampling.from(), inSampling.to(), inSampling.step(), inSampling.stride());
  Sampling1D<typename TRasterOut::Value> xConvolvedSampling(xConvolvedWidth, xConvolved.data());
  printf(
      "xConvolvedSampling: %i-%i:%i:%i\n",
      xConvolvedSampling.from(),
      xConvolvedSampling.to(),
      xConvolvedSampling.step(),
      xConvolvedSampling.stride());
  for (long y = yFrom; y <= yTo; ++y) {
    printf(
        "in: %i, xConvolved: %i\n",
        std::distance(in.data(), inSampling.data()),
        std::distance(xConvolved.data(), xConvolvedSampling.data()));
    convolve1DZero(inSampling, kernel, xConvolvedSampling);
    inSampling.data(inSampling.data() + in.shape()[0]);
    xConvolvedSampling.data(xConvolvedSampling.data() + xConvolvedWidth);
    for (long x = 0; x < xConvolvedWidth; ++x)
      printf("%i, %i: %i\n", x, y, xConvolved[{x, y - yFrom}]);
  }

  printf("xConvolved:\n%ix%i\n", xConvolvedWidth, xConvolvedHeight);
  for (long y = 0; y < xConvolvedHeight; ++y) {
    for (long x = 0; x < xConvolvedWidth; ++x) {
      printf("%i ", xConvolved[{x, y}]);
    }
    printf("\n");
  }

  // Convolve along y-axis
  Sampling1D<typename TRasterOut::Value> ySampling(xConvolved.template length<1>(), xConvolved.data());
  ySampling.from(region.front[1]).to(region.back[1]).step(yStep).stride(xConvolvedWidth);
  printf("ySampling: %i-%i:%i:%i\n", ySampling.from(), ySampling.to(), ySampling.step(), ySampling.stride());
  Sampling1D<typename TRasterOut::Value> outSampling(out.template length<1>(), out.data());
  outSampling.stride(out.template length<0>());
  printf("outSampling: %i-%i:%i:%i\n", outSampling.from(), outSampling.to(), outSampling.step(), outSampling.stride());
  for (long x = xFrom; x <= xTo; x += xStep) {
    printf(
        "xConvolved: %i, out: %i\n",
        std::distance(xConvolved.data(), ySampling.data()),
        std::distance(out.data(), outSampling.data()));
    convolve1DZero(ySampling, kernel, outSampling);
    ySampling.data(ySampling.data() + 1);
    outSampling.data(outSampling.data() + 1);
    for (long y = 0; y < out.template length<1>(); ++y)
      printf("%i, %i: %i\n", x, y, out[{x, y}]);
  }
}

template <typename TRasterIn, typename TKernel, typename TRasterOut>
void convolveXyZero(const TRasterIn& in, const TKernel& kernel, TRasterOut& out) {
  convolveXyZero(in, in.domain(), Position<2>::one(), kernel, out);
}

} // namespace Image2D
} // namespace Phi

#endif // _PHIBOX_IMAGEPROCESSING_H
