// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef _PHIDUFFIEUX_MONOCHROMATICOPTICS_H
#define _PHIDUFFIEUX_MONOCHROMATICOPTICS_H

#include "PhiFourier/Dft.h"

#include <vector>

namespace Phi {
namespace Duffieux {

/**
 * @brief Monochromatic data buffers and transforms.
 */
class MonochromaticOptics {

  friend class MonochromaticSystem;

public:
  /**
   * @brief Constructor.
   * @param lambda The wavelength
   * @param diameter The pupil mask diameter
   * @param alphas The Zernike coefficients
   */
  MonochromaticOptics(double lambda, long diameter, std::vector<double> alphas) :
      m_minusTwoPiOverLambda(-2 * 3.1415926 / lambda), m_alphas(std::move(alphas)), m_pupilToPsf({diameter, diameter}),
      m_pupilAmplitude(m_pupilToPsf.in()), m_psfAmplitude(m_pupilToPsf.out()), m_psfIntensity({diameter, diameter}) {}

  /**
   * @brief Get the pupil amplitude.
   */
  const Fourier::ComplexDftBuffer& pupilAmplitude() const {
    return m_pupilAmplitude;
  }

  /**
   * @brief Get the PSF amplitude.
   */
  const Fourier::ComplexDftBuffer& psfAmplitude() const {
    return m_psfAmplitude;
  }

  /**
   * @brief Get the PSF intensity.
   */
  const Fourier::RealDftBuffer& psfIntensity() const {
    return m_psfIntensity;
  }

  /**
   * @brief Evaluate the pupil amplitude from the pupil mask and Zernike coefficients.
   * @param mask The pupil mask raster
   * @param zernikes The Zernike basis ordered as (i, x, y)
   */
  template <typename TMask, typename TZernikes>
  Fourier::ComplexDftBuffer& evalPupilAmplitude(const TMask& mask, const TZernikes& zernikes) {
    auto maskData = mask.data();
    auto zernikesData = zernikes.data();
    const auto size = m_alphas.size();
    for (auto it = m_pupilAmplitude.begin(); it != m_pupilAmplitude.end(); ++it, ++maskData, zernikesData += size) {
      if (*maskData != 0) {
        *it = evalPhase(*maskData, zernikesData);
      } else {
        *it = 0;
      }
    }
    return m_pupilAmplitude;
  }

  /**
   * @brief Evaluate the PSF amplitude as the DFT of the pupil amplitude.
   */
  Fourier::ComplexDftBuffer& evalPsfAmplitude() {
    m_pupilToPsf.transform();
    return m_psfAmplitude;
  }

  /**
   * @brief Evaluate the PSF intensity as the norm squared of the PSF amplitude.
   */
  Fourier::RealDftBuffer& evalPsfIntensity() {
    m_psfIntensity.generate(
        [](const auto& amp) {
          return std::norm(amp);
        },
        m_psfAmplitude);
    return m_psfIntensity;
  }

private:
  /**
   * @brief Evaluate some local phase.
   * @param mask The local pupil mask value
   * @param zernikes The local Zernike basis values
   */
  std::complex<double> evalPhase(double mask, const double* zernikes) {
    double sum = 0;
    auto zIt = zernikes;
    for (auto aIt = m_alphas.begin(); aIt != m_alphas.end(); ++aIt, ++zIt) {
      sum += (*aIt) * (*zIt);
    }
    return mask * std::exp(std::complex<double>(0, m_minusTwoPiOverLambda * sum));
  }

  double m_minusTwoPiOverLambda;
  std::vector<double> m_alphas;
  Fourier::ComplexDft m_pupilToPsf;
  Fourier::ComplexDftBuffer& m_pupilAmplitude;
  Fourier::ComplexDftBuffer& m_psfAmplitude;
  Fourier::RealDftBuffer m_psfIntensity;
};

} // namespace Duffieux
} // namespace Phi

#endif // _PHIDUFFIEUX_MONOCHROMATICOPTICS_H
