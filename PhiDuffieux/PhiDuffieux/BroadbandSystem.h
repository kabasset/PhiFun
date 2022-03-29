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

struct SystemTfList {
  using Prerequisite = void;
  using Return = std::list<const Fourier::ComplexDftBuffer*>;
};

struct BroadbandTf {
  using Prerequisite = SystemTfList;
  using Return = const Fourier::ComplexDftBuffer&;
};

struct BroadbandPsf {
  using Prerequisite = BroadbandTf;
  using Return = const Fourier::RealDftBuffer&;
};

class BroadbandSystem : public Framework::StepperPipeline<BroadbandSystem> {

public:
  struct Params {
    std::vector<double> wavelengths;
    std::vector<double> integrationWavelengths;
    std::vector<double> spectrum;
    Fourier::Position shape;
  };

  BroadbandSystem(MonochromaticOptics::Params optics, MonochromaticSystem::Params system, Params params) :
      m_params(std::move(params)), m_systems(), m_tfs(),
      m_integrator(m_params.wavelengths, m_params.integrationWavelengths), m_tfToPsf(m_params.shape) {
    m_systems.reserve(m_params.wavelengths.size());
    for (auto lambda : m_params.wavelengths) {
      auto mop = optics;
      mop.updateWavelength(lambda);
      m_systems.emplace_back(std::move(mop), system);
    }
  }

protected:
  template <typename S>
  typename S::Return doGet();

  template <typename S>
  void doEvaluate();

private:
  Params m_params;
  std::vector<MonochromaticSystem> m_systems;
  typename SystemTfList::Return m_tfs;
  Spline::SplineIntegrator m_integrator;
  Fourier::RealDft::Inverse m_tfToPsf;
};

template <>
typename SystemTfList::Return BroadbandSystem::doGet<SystemTfList>() {
  return m_tfs;
}

template <>
void BroadbandSystem::doEvaluate<SystemTfList>() {
  for (auto& s : m_systems) {
    m_tfs.push_back(&s.get<SystemTf>());
  }
}

template <>
typename BroadbandTf::Return BroadbandSystem::doGet<BroadbandTf>() {
  return m_tfToPsf.in();
}

template <>
void BroadbandSystem::doEvaluate<BroadbandTf>() {
  const auto& inShape = m_systems[0].optics().params().shape;
  const auto& outShape = m_systems[0].params().shape;
  const auto width = outShape[0];
  const auto height = outShape[1];
  const double xFactor = double(inShape[0] - 1) / (width - 1);
  const double yFactor = double(inShape[1] - 1) / (height - 1);
  const auto depth = m_systems.size();
  std::vector<std::complex<double>> y(depth);
  std::vector<std::complex<double>> z(depth);
  auto* it = m_tfToPsf.in().begin();
  for (long j = 0; j < height; ++j) {
    const double v = yFactor * j;
    for (long i = 0; i < width; ++i, ++it) {
      const double u = xFactor * i;
      for (std::size_t l = 0; l < depth; ++l) {
        const auto lambda = m_params.wavelengths[l];
        y[l] = Image2D::bilinear<std::complex<double>>(m_systems[l].get<SystemTf>(), u / lambda, v / lambda);
      }
      z = m_integrator.knotZ(y.data());
      *it = m_integrator.integrate(y.data(), z.data(), m_params.spectrum.data());
    }
  }
}

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
