// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef _PHIDUFFIEUX_MONOCHROMATICSYSTEM_H
#define _PHIDUFFIEUX_MONOCHROMATICSYSTEM_H

#include "PhiBox/ImageProcessing.h"
#include "PhiDuffieux/MonochromaticOptics.h"

namespace Phi {
namespace Duffieux {

/**
 * @brief Monochromatic system model.
 */
class MonochromaticSystem : public Framework::StepperPipeline<MonochromaticSystem> {

public:
  /**
   * @brief Distortion model.
   */
  struct Distortion {
    double ux, uy, vx, vy; ///< The distortion parameters
  };

  /**
   * @brief Non optical parameters.
   */
  struct Parameters {
    Fourier::Position shape; ///< The broadband logical data shape
    const std::complex<double>* nonOpticalTf; ///< The non optical transfer function
    Distortion distortion; ///< The distortion model
  };

  /**
   * @brief Constructor.
   */
  MonochromaticSystem(MonochromaticOptics::Parameters opticalParameters, Parameters nonOpticalParameters) :
      m_parameters(std::move(nonOpticalParameters)), m_optics(std::move(opticalParameters)),
      m_psfToTf(m_optics.m_parameters.shape, m_optics.m_psfIntensity.data()), m_tfToPsf(m_parameters.shape) {}

  /**
   * @brief Get the parameters.
   */
  const Parameters& parameters() const {
    return m_parameters;
  }

  /**
   * @brief Get the optical model.
   */
  const MonochromaticOptics& optics() {
    return m_optics;
  }

  /**
   * @brief Get the wavelength.
   */
  double wavelength() const {
    return m_optics.wavelength();
  }

  /**
   * @brief Update the wavelength and optionally the associated parameters.
   * @param wavelength The new wavelength
   * @param coefficients The new Zernike coefficients
   * @param nonOpticalTf The new non optical TF
   * @warning
   * This method resets the pipeline.
   */
  void
  update(double wavelength, const double* coefficients = nullptr, const std::complex<double>* nonOpticalTf = nullptr) {
    m_optics.update(wavelength, coefficients);
    if (nonOpticalTf) {
      m_parameters.nonOpticalTf = nonOpticalTf;
    }
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
  Parameters m_parameters; ///< The system parameters
  MonochromaticOptics m_optics; ///< The optical model
  Fourier::RealDft m_psfToTf; ///< The unwarped PSF to TF transform
  Fourier::RealDft::Inverse m_tfToPsf; ///< The warped TF to PSF transform
};

/**
 * @brief The system TF without distortion in the monochromatic grid.
 */
struct SystemTf : Framework::PipelineStep<PsfIntensity, const Fourier::ComplexDftBuffer&> {};

template <>
inline const Fourier::ComplexDftBuffer& MonochromaticSystem::doGet<SystemTf>() {
  return m_psfToTf.out();
}

template <>
inline void MonochromaticSystem::doEvaluate<SystemTf>() {
  m_psfToTf.transform();
  const Euclid::Fits::PtrRaster<const std::complex<double>> noTf(m_psfToTf.outShape(), m_parameters.nonOpticalTf);
  m_psfToTf.out() *= noTf;
}

/**
 * @brief The regrided and distorted TF.
 */
struct WarpedSystemTf : Framework::PipelineStep<SystemTf, const Fourier::ComplexDftBuffer&> {};

template <>
inline const Fourier::ComplexDftBuffer& MonochromaticSystem::doGet<WarpedSystemTf>() {
  return m_tfToPsf.in();
}

template <>
inline void MonochromaticSystem::doEvaluate<WarpedSystemTf>() {
  const auto width = m_tfToPsf.in().shape()[0];
  const auto height = m_tfToPsf.in().shape()[1];
  const auto wavelength = m_optics.m_parameters.wavelength;
  const double xFactor = wavelength * (m_parameters.shape[0] - 1) / (width - 1);
  const double yFactor = wavelength * (m_parameters.shape[1] - 1) / (height - 1);
  const double uxFactor = m_parameters.distortion.ux * xFactor;
  const double uyFactor = m_parameters.distortion.uy * yFactor;
  const double vxFactor = m_parameters.distortion.vx * xFactor;
  const double vyFactor = m_parameters.distortion.vy * yFactor;
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

/**
 * @brief The regrided and distorted PSF.
 */
struct WarpedSystemPsf : Framework::PipelineStep<WarpedSystemTf, const Fourier::RealDftBuffer&> {};

template <>
inline const Fourier::RealDftBuffer& MonochromaticSystem::doGet<WarpedSystemPsf>() {
  return m_tfToPsf.out();
}

template <>
inline void MonochromaticSystem::doEvaluate<WarpedSystemPsf>() {
  m_tfToPsf.transform().normalize();
}

} // namespace Duffieux
} // namespace Phi

#endif // _PHIDUFFIEUX_MONOCHROMATICSYSTEM_H
