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

} // namespace Fourier
} // namespace Phi

#endif // _PHIFOURIER_DFT_H
