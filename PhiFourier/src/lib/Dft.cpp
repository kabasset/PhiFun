// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "PhiFourier/Dft.h"

#include "PhiFourier/DftMemory.h"

namespace Phi {
namespace Fourier {

template <>
Position RealDftType::Parent::outShape(const Position& shape) {
  return {shape[0] / 2 + 1, shape[1]};
}

template <>
fftw_plan Internal::allocateFftwPlan<RealDftType>(DftBuffer<double>& in, DftBuffer<std::complex<double>>& out) {
  const auto& shape = in.shape();
  const int width = static_cast<int>(shape[0]);
  const int height = static_cast<int>(shape[1]);
  int n[] = {height, width}; // FFTW ordering
  return fftw_plan_many_dft_r2c(
      2, // rank
      n, // n
      shape[2], // howmany
      reinterpret_cast<double*>(in.data()), // in // FIXME reinterpret needed?
      nullptr, // inembed
      1, // istride
      width * height, // idist
      reinterpret_cast<fftw_complex*>(out.data()), // out
      nullptr, // onembed
      1, // ostride
      (width / 2 + 1) * height, // odist
      FFTW_MEASURE); // FIXME other flags?
}

template <>
fftw_plan
Internal::allocateFftwPlan<Inverse<RealDftType>>(DftBuffer<std::complex<double>>& in, DftBuffer<double>& out) {
  const auto& shape = out.shape();
  const int width = static_cast<int>(shape[0]);
  const int height = static_cast<int>(shape[1]);
  int n[] = {height, width}; // FFTW ordering
  return fftw_plan_many_dft_c2r(
      2, // rank
      n, // n
      shape[2], // howmany
      reinterpret_cast<fftw_complex*>(in.data()), // in
      nullptr, // inembed
      1, // istride
      (width / 2 + 1) * height, // idist
      reinterpret_cast<double*>(out.data()), // out // FIXME reinterpret needed?
      nullptr, // onembed
      1, // ostride
      width * height, // odist
      FFTW_MEASURE); // FIXME other flags?
}

template <>
fftw_plan
Internal::allocateFftwPlan<ComplexDftType>(DftBuffer<std::complex<double>>& in, DftBuffer<std::complex<double>>& out) {
  const auto& shape = in.shape();
  const int width = static_cast<int>(shape[0]);
  const int height = static_cast<int>(shape[1]);
  int n[] = {height, width}; // FFTW ordering
  return fftw_plan_many_dft(
      2, // rank
      n,
      shape[2], // howmany,
      reinterpret_cast<fftw_complex*>(in.data()), // in
      nullptr, // inembed
      1, // istride
      width * height, // idist
      reinterpret_cast<fftw_complex*>(out.data()), // out
      nullptr, // onembed
      1, // ostride
      width * height, // odist
      FFTW_FORWARD, // sign
      FFTW_MEASURE); // FIXME other flags?
}

template <>
fftw_plan Internal::allocateFftwPlan<Inverse<ComplexDftType>>(
    DftBuffer<std::complex<double>>& in,
    DftBuffer<std::complex<double>>& out) {
  const auto& shape = out.shape();
  const int width = static_cast<int>(shape[0]);
  const int height = static_cast<int>(shape[1]);
  int n[] = {height, width}; // FFTW ordering
  return fftw_plan_many_dft(
      2, // rank
      n,
      shape[2], // howmany,
      reinterpret_cast<fftw_complex*>(in.data()), // in
      nullptr, // inembed
      1, // istride
      width * height, // idist
      reinterpret_cast<fftw_complex*>(out.data()), // out
      nullptr, // onembed
      1, // ostride
      width * height, // odist
      FFTW_BACKWARD, // sign
      FFTW_MEASURE); // FIXME other flags?
}

template <>
Position HermitianComplexDftType::Parent::inShape(const Position& shape) {
  return {shape[0] / 2 + 1, shape[1]};
}

template <>
Position HermitianComplexDftType::Parent::outShape(const Position& shape) {
  return {shape[0] / 2 + 1, shape[1]};
}

template <>
fftw_plan Internal::allocateFftwPlan<HermitianComplexDftType>(
    DftBuffer<std::complex<double>>& in,
    DftBuffer<std::complex<double>>& out) {
  return Internal::allocateFftwPlan<ComplexDftType>(in, out);
}

template <>
fftw_plan Internal::allocateFftwPlan<Inverse<HermitianComplexDftType>>(
    DftBuffer<std::complex<double>>& in,
    DftBuffer<std::complex<double>>& out) {
  return Internal::allocateFftwPlan<Inverse<ComplexDftType>>(in, out);
}

} // namespace Fourier
} // namespace Phi