// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef _PHIDUFFIEUX_MONOCHROMATICSYSTEM_H
#define _PHIDUFFIEUX_MONOCHROMATICSYSTEM_H

#include "PhiBox/ImageProcessing.h"
#include "PhiDuffieux/MonochromaticOptics.h"

namespace Phi {
namespace Duffieux {

struct SystemTf {
  using Prerequisite = PsfIntensity;
  using Return = const Fourier::ComplexDftBuffer&;
};

struct WarpedSystemTf {
  using Prerequisite = SystemTf;
  using Return = const Fourier::ComplexDftBuffer&;
};

struct WarpedSystemPsf {
  using Prerequisite = WarpedSystemTf;
  using Return = const Fourier::RealDftBuffer&;
};

/**
 * @brief Monochromatic system model.
 */
class MonochromaticSystem : public Framework::StepperPipeline<MonochromaticSystem> {

public:
  /**
   * @brief Non optical parameters.
   */
  struct Params {
    Fourier::Position shape; ///< The broadband logical data shape
    const Fourier::ComplexDftBuffer& nonOpticalTf; ///< The non optical transfer function
    double distortion[4]; ///< The distortion coefficients ordered as: `ux, ux, vx, vy`
  };

  /**
   * @brief Constructor.
   */
  MonochromaticSystem(MonochromaticOptics::Params opticalParams, Params nonOpticalParams) :
      m_params(std::move(nonOpticalParams)), m_optics(std::move(opticalParams)),
      m_psfToTf(m_optics.m_params.shape, m_optics.m_psfIntensity.data()), m_tfToPsf(m_params.shape) {}

  /**
   * @brief Get the optical model.
   */
  const MonochromaticOptics& optics() {
    return m_optics;
  }

  void updateLambda(double lambda) {
    m_optics.updateLambda(lambda);
    reset();
  }

protected:
  /**
   * @brief Get the result of step `S`.
   */
  template <typename S>
  typename S::Return doGet() {
    return m_optics.template get<S>();
  }

  /**
   * @brief Evaluate step `S`.
   */
  template <typename S>
  void doEvaluate() {
    m_optics.template get<S>();
  }

private:
  Params m_params; ///< The system parameters
  MonochromaticOptics m_optics; ///< The optical model
  Fourier::RealDft m_psfToTf; ///< The unwarped PSF to TF transform
  Fourier::RealDft::Inverse m_tfToPsf; ///< The warped TF to PSF transform
};

template <>
inline const Fourier::ComplexDftBuffer& MonochromaticSystem::doGet<SystemTf>() {
  return m_psfToTf.out();
}

template <>
inline const Fourier::ComplexDftBuffer& MonochromaticSystem::doGet<WarpedSystemTf>() {
  return m_tfToPsf.in();
}

template <>
inline const Fourier::RealDftBuffer& MonochromaticSystem::doGet<WarpedSystemPsf>() {
  return m_tfToPsf.out();
}

template <>
inline void MonochromaticSystem::doEvaluate<SystemTf>() {
  m_psfToTf.transform();
  m_psfToTf.out() *= m_params.nonOpticalTf;
}

template <>
inline void MonochromaticSystem::doEvaluate<WarpedSystemTf>() {
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
inline void MonochromaticSystem::doEvaluate<WarpedSystemPsf>() {
  m_tfToPsf.transform().normalize();
}

} // namespace Duffieux
} // namespace Phi

#endif // _PHIDUFFIEUX_MONOCHROMATICSYSTEM_H
