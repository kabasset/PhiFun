// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef _PHIDUFFIEUX_BROADBANDSYSTEM_H
#define _PHIDUFFIEUX_BROADBANDSYSTEM_H

#include "PhiBox/SplineIntegrator.h"
#include "PhiBox/StepperPipeline.h"
#include "PhiDuffieux/MonochromaticSystem.h"
#include "PhiFourier/Dft.h"

#include <list>

namespace Phi {
namespace Duffieux {

class BroadbandSystem : public Framework::StepperPipeline<BroadbandSystem> {

public:
  using Stack = Euclid::Fits::VecRaster<std::complex<double>, 3>;
  struct Params {
    std::vector<double> wavelengths;
    std::vector<double> integrationWavelengths;
    std::vector<double> spectrum;
    Fourier::Position shape;
  };

  BroadbandSystem(MonochromaticOptics::Params optics, MonochromaticSystem::Params system, Params params) :
      m_params(std::move(params)), // Broadband parameters
      m_system(std::move(optics), std::move(system)), // Monochromatic parameters
      m_tfs(
          {m_system.optics().params().shape[0] / 2 + 1,
           m_system.optics().params().shape[1],
           long(m_params.wavelengths.size())}), // TF stack
      m_integrator(m_params.wavelengths, m_params.integrationWavelengths), // Spline integrator
      m_tfToPsf(m_params.shape) {} // iDFT

protected:
  template <typename S>
  typename S::Return doGet();

  template <typename S>
  void doEvaluate();

private:
  Params m_params;
  MonochromaticSystem m_system;
  Stack m_tfs;
  Spline::SplineIntegrator m_integrator;
  Fourier::RealDft::Inverse m_tfToPsf;
};

/**
 * @brief The stack of monochromatic system TFs.
 */
struct SystemTfStack : Framework::PipelineStep<void, const BroadbandSystem::Stack&> {};

template <>
typename SystemTfStack::Return BroadbandSystem::doGet<SystemTfStack>() {
  return m_tfs;
}

template <>
void BroadbandSystem::doEvaluate<SystemTfStack>() {
  auto* it = m_tfs.data();
  for (auto& lambda : m_params.wavelengths) {
    m_system.updateWavelength(lambda);
    const auto& tf = m_system.get<SystemTf>();
    std::copy(tf.begin(), tf.end(), it);
    it += shapeSize(m_params.shape);
  }
}

/**
 * @brief The broadband TF.
 */
struct BroadbandTf : Framework::PipelineStep<SystemTfStack, const Fourier::ComplexDftBuffer&> {};

template <>
typename BroadbandTf::Return BroadbandSystem::doGet<BroadbandTf>() {
  return m_tfToPsf.in();
}

template <>
void BroadbandSystem::doEvaluate<BroadbandTf>() {
  const auto& inShape = m_tfs.shape();
  const auto& outShape = m_tfToPsf.inShape();
  const auto width = outShape[0];
  const auto height = outShape[1];
  const double xFactor = double(inShape[0] - 1) / (width - 1);
  const double yFactor = double(inShape[1] - 1) / (height - 1);
  const auto depth = inShape[2];
  std::vector<std::complex<double>> y(depth);
  std::vector<std::complex<double>> z(depth);
  auto* it = m_tfToPsf.in().begin();
  for (long j = 0; j < height; ++j) {
    const double v = yFactor * j;
    for (long i = 0; i < width; ++i, ++it) {
      const double u = xFactor * i;
      for (long l = 0; l < depth; ++l) {
        const auto lambda = m_params.wavelengths[l];
        y[l] = Image2D::bilinear<std::complex<double>>(m_tfs, u / lambda, v / lambda, l);
      }
      z = m_integrator.knotZ(y.data());
      *it = m_integrator.integrate(y.data(), z.data(), m_params.spectrum.data());
    }
  }
}

/**
 * @brief The broadband PSF.
 */
struct BroadbandPsf : Framework::PipelineStep<BroadbandTf, const Fourier::RealDftBuffer&> {};

template <>
typename BroadbandPsf::Return BroadbandSystem::doGet<BroadbandPsf>() {
  return m_tfToPsf.out();
}

template <>
void BroadbandSystem::doEvaluate<BroadbandPsf>() {
  m_tfToPsf.transform().normalize();
}

} // namespace Duffieux
} // namespace Phi

#endif // _PHIDUFFIEUX_BROADBANDSYSTEM_H
