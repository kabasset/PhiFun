// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "PhiDuffieux/MonochromaticSystem.h"

namespace Phi {
namespace Duffieux {

template <>
const Fourier::ComplexDftBuffer& MonochromaticSystem::doGet<SystemTf>() const {
  return m_psfToTf.out();
}

template <>
const Fourier::ComplexDftBuffer& MonochromaticSystem::doGet<WarpedSystemTf>() const {
  return m_tfToPsf.in();
}

template <>
const Fourier::RealDftBuffer& MonochromaticSystem::doGet<WarpedSystemPsf>() const {
  return m_tfToPsf.out();
}

template <>
void MonochromaticSystem::doEvaluate<SystemTf>() {
  m_psfToTf.transform();
  m_psfToTf.out() *= m_params.nonOpticalTf;
}

template <>
void MonochromaticSystem::doEvaluate<WarpedSystemTf>() {
  const auto width = m_tfToPsf.in().shape()[0];
  const auto height = m_tfToPsf.in().shape()[1];
  const auto lambda = m_optics.m_params.wavelength;
  const double xFactor = lambda * (m_params.shape[0] - 1) / (width - 1);
  const double yFactor = lambda * (m_params.shape[1] - 1) / (height - 1);
  const double uxFactor = m_params.distortion[0] * xFactor;
  const double uyFactor = m_params.distortion[1] * yFactor;
  const double vxFactor = m_params.distortion[2] * xFactor;
  const double vyFactor = m_params.distortion[3] * yFactor;
  auto* it = m_tfToPsf.in().begin();
  auto& stf = m_psfToTf.out();
  for (long y = 0; y < height; ++y) {
    const double uyValue = uyFactor * y;
    const double vyValue = vyFactor * y;
    for (long x = 0; x < width; ++x, ++it) {
      const double u = uxFactor * x + uyValue;
      const double v = vxFactor * x + vyValue;
      *it = Image2D::bilinear<std::complex<double>>(stf, u, v);
    }
  }
}

template <>
void MonochromaticSystem::doEvaluate<WarpedSystemPsf>() {
  m_tfToPsf.transform().normalize();
}

} // namespace Duffieux
} // namespace Phi
