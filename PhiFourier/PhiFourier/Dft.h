// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef _PHIFOURIER_DFT_H
#define _PHIFOURIER_DFT_H

#include "PhiFourier/DftPlan.h"

namespace Phi {
namespace Fourier {

/**
 * @brief Inverse type of a `DftType`.
 */
template <typename TType>
struct Inverse; // Forward declaration for DftType

/**
 * @brief Base DFT type to be inherited.
 */
template <typename TType, typename TIn, typename TOut>
struct DftType {

  /**
   * @brief The parent `DftType` class.
   */
  using Parent = DftType;

  /**
   * @brief The type tag.
   */
  using Type = TType;

  /**
   * @brief The input value type.
   */
  using InValue = TIn;

  /**
   * @brief The output value type.
   */
  using OutValue = TOut;

  /**
   * @brief The tag of the inverse transform type.
   */
  using InverseType = Inverse<TType>;

  /**
   * @brief Input buffer shape.
   * @param shape The logical shape
   */
  static Position inShape(const Position& shape) {
    return shape;
  }

  /**
   * @brief Output buffer shape.
   * @param shape The logical shape
   */
  static Position outShape(const Position& shape) {
    return shape;
  }
};

/**
 * @brief Specialization for inverse types.
 */
template <typename TType, typename TIn, typename TOut>
struct DftType<Inverse<TType>, TIn, TOut> {

  using Parent = DftType;
  using Type = Inverse<TType>;
  using InValue = TOut;
  using OutValue = TIn;
  using InverseType = TType;

  static Position inShape(const Position& shape) {
    return DftType<TType, TIn, TOut>::outShape(shape);
  }

  static Position outShape(const Position& shape) {
    return DftType<TType, TIn, TOut>::inShape(shape);
  }
};

template <typename TType>
struct Inverse : DftType<Inverse<TType>, typename TType::InValue, typename TType::OutValue> {};

/**
 * @brief Real DFT type.
 */
struct RealDftType;
struct RealDftType : DftType<RealDftType, double, std::complex<double>> {};
template <>
Position RealDftType::Parent::outShape(const Position& shape);

/**
 * @brief Complex DFT type.
 */
struct ComplexDftType;
struct ComplexDftType : DftType<ComplexDftType, std::complex<double>, std::complex<double>> {};

/**
 * @brief Complex DFT type with Hermitian symmertry.
 */
struct HermitianComplexDftType;
struct HermitianComplexDftType : DftType<HermitianComplexDftType, std::complex<double>, std::complex<double>> {};
template <>
Position HermitianComplexDftType::Parent::inShape(const Position& shape);
template <>
Position HermitianComplexDftType::Parent::outShape(const Position& shape);

/**
 * @brief Real DFT plan.
 */
using RealDft = DftPlan<RealDftType>;

/**
 * @brief Complex DFT plan.
 */
using ComplexDft = DftPlan<ComplexDftType>;

/**
 * @brief Complex DFT plan with Hermitian symmetry.
 */
using HermitianComplexDft = DftPlan<HermitianComplexDftType>;

/// @cond
namespace Internal {

template <typename TIter>
void swapRanges(TIter aBegin, TIter aEnd, TIter bBegin) {
  TIter aIt = aBegin;
  TIter bIt = bBegin;
  while (aIt != aEnd) {
    std::iter_swap(aIt++, bIt++);
  }
}

} // namespace Internal
/// @endcond

/**
 * @brief Swap the quadrants of a raster in place.
 */
template <typename TRaster>
TRaster& shiftDft(TRaster& raster) {
  const auto width = raster.shape()[0];
  const auto height = raster.shape()[1];
  if (width % 2 != 0 || height % 2 != 0) {
    throw std::runtime_error("shiftDft() only works with even sizes as of today.");
  }
  const long halfWidth = raster.shape()[0] / 2;
  const long halfHeight = raster.shape()[1] / 2;

  for (long y = 0; y < halfHeight; ++y) {

    // Swap UL with LR
    auto ulBegin = &raster[{0, y}];
    auto ulEnd = ulBegin + halfWidth;
    auto lrBegin = &raster[{halfWidth, y + halfHeight}];
    Internal::swapRanges(ulBegin, ulEnd, lrBegin);

    // Swap UR and LL
    auto urBegin = ulEnd;
    auto urEnd = urBegin + halfWidth;
    auto llBegin = &raster[{0, y + halfHeight}];
    Internal::swapRanges(urBegin, urEnd, llBegin);
  }

  return raster;
}

/**
 * @brief Apply the norm squared to each element of a raster.
 */
template <typename TComplexRaster, typename TRealRaster>
TRealRaster& norm2(const TComplexRaster& input, TRealRaster& output) {
  output.generate(
      [](const std::complex<double>& c) {
        return std::norm(c);
      },
      input);
  return output;
}

/**
 * @copybrief norm2()
 */
template <typename TComplexRaster>
Euclid::Fits::VecRaster<double> norm2(const TComplexRaster& input) {
  Euclid::Fits::VecRaster<double> output(input.shape());
  norm2(input, output);
  return output;
}

} // namespace Fourier
} // namespace Phi

#endif // _PHIFOURIER_DFT_H
