// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef _PHIBOX_CONVOLUTION_H
#define _PHIBOX_CONVOLUTION_H

namespace Phi {

template <typename T, long N, typename TBoundary>
class SeparableKernel {
public:
  using Value = T;

  template <typename TKernel>
  SeparableKernel(const TKernel& coefficients, const Position<N>& center, TBoundary boundary = {});

  template <typename TRasterIn>
  VecRaster<typename TRasterIn::Value, TRasterIn::Dim> operator*(const TRasterIn& in) const;

  template <typename T, typename TRasterIn>
  VecRaster<T, TRasterIn::Dim> as(const TRasterIn& in) const;

  template <typename TRasterIn, typename TRasterOut>
  void to(const TRasterIn& in, TRasterOut& out) const;

  template <typename T, typename TRasterIn>
  VecRaster<T, TRasterIn::Dim> sparseAs(const TRasterIn& in, const LinearSampling<Dim>& sampling) const;

  template <typename TRasterIn, typename TRasterOut>
  void sparseTo(const TRasterIn& in, const LinearSampling<Dim>& sampling, TRasterOut& out) const;

private:
  void processZero(const TRasterIn& in, const Region<Dim>& region, const Position<Dim>& step, TRasterOut& out) const;

  template <typename UBoundary = TBoundary, std::enable_if_t<...>* = nullptr>
  void postProcessBoundaries() const;

  VecRaster<T, N> m_coefficients;
  Position<N> m_center;
  TBoundary m_boundary;
};

template <typename TRaster, typename T, long N, typename TBoundary>
VecRaster<typename TRaster::Value, TRaster::N>
operator*(const TRaster& in, const SeparableKernel<T, N, TBoundary>& kernel) {
  return kernel * in;
}

} // namespace Phi

#endif // _PHIBOX_CONVOLUTION_H
