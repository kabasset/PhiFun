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
 * @brief Monochromatic data buffers and transforms with lazy evaluation.
 */
class MonochromaticOptics {
  friend class MonochromaticSystem;

  /**
   * @brief 3.14 ;)
   */
  static constexpr double m_pi = std::acos(-1.);

public:
  /**
   * @brief Constructor.
   * @param lambda The wavelength
   * @param diameter The pupil mask
   * @param basis The Zernike basis
   * @param alphas The Zernike coefficients
   */
  template <typename TMask, typename TZernikes>
  MonochromaticOptics(double lambda, const TMask& mask, const TZernikes& basis, std::vector<double> alphas) :
      m_state(Initial), m_wavenumber(2 * m_pi / lambda), m_maskData(mask.data()), m_zernikesData(basis.data()),
      m_alphas(std::move(alphas)), m_pupilToPsf(mask.shape()), m_pupilAmplitude(m_pupilToPsf.in()),
      m_psfAmplitude(m_pupilToPsf.out()), m_psfIntensity(mask.shape()) {}

  /**
   * @brief Update the wavelength.
   */
  void updateLambda(double lambda) {
    m_wavenumber = 2 * m_pi / lambda;
    m_state = State::Initial;
  }

  /**
   * @brief Get the pupil amplitude.
   */
  const Fourier::ComplexDftBuffer& pupilAmplitude() {
    if (m_state < State::PupilAmplitude) {
      evalPupilAmplitude();
    }
    return m_pupilAmplitude;
  }

  /**
   * @brief Get the PSF amplitude.
   */
  const Fourier::ComplexDftBuffer& psfAmplitude() {
    if (m_state < State::PsfAmplitude) {
      evalPsfAmplitude();
    }
    return m_psfAmplitude;
  }

  /**
   * @brief Get the PSF intensity.
   */
  const Fourier::RealDftBuffer& psfIntensity() {
    if (m_state < State::PsfIntensity) {
      evalPsfIntensity();
    }
    return m_psfIntensity;
  }

protected:
  /**
   * @brief The current computation state.
   */
  enum State
  {
    Initial,
    PupilAmplitude,
    PsfAmplitude,
    PsfIntensity
  };

  /**
   * @brief Evaluate the pupil amplitude from the pupil mask and Zernike coefficients.
   * @param mask The pupil mask raster
   * @param zernikes The Zernike basis ordered as (i, x, y)
   */
  Fourier::ComplexDftBuffer& evalPupilAmplitude() {
    const auto size = m_alphas.size();
    const double* maskIt = m_maskData;
    const double* zernikesIt = m_zernikesData;
    for (auto it = m_pupilAmplitude.begin(); it != m_pupilAmplitude.end(); ++it, ++maskIt, zernikesIt += size) {
      if (*maskIt != 0) {
        *it = evalPhase(*maskIt, zernikesIt);
      } else {
        *it = 0;
      }
    }
    m_state = State::PupilAmplitude;
    return m_pupilAmplitude;
  }

  /**
   * @brief Evaluate the PSF amplitude as the DFT of the pupil amplitude.
   */
  Fourier::ComplexDftBuffer& evalPsfAmplitude() {
    if (m_state < State::PupilAmplitude) {
      evalPupilAmplitude();
    }
    m_pupilToPsf.transform().normalize();
    m_state = State::PsfAmplitude;
    return m_psfAmplitude;
  }

  /**
   * @brief Evaluate the PSF intensity as the norm squared of the PSF amplitude.
   */
  Fourier::RealDftBuffer& evalPsfIntensity() {
    if (m_state < State::PsfAmplitude) {
      evalPsfAmplitude();
    }
    norm2(m_psfAmplitude, m_psfIntensity);
    m_state = State::PsfIntensity;
    return m_psfIntensity;
  }

private:
  /**
   * @brief Evaluate some local phase.
   * @param mask The local pupil mask value
   * @param zernikes The local Zernike basis values
   */
  std::complex<double> evalPhase(double mask, const double* zernikes) {
    double minusPhi = 0;
    auto zIt = zernikes;
    for (auto aIt = m_alphas.begin(); aIt != m_alphas.end(); ++aIt, ++zIt) {
      minusPhi -= (*aIt) * (*zIt);
    }
    return mask * std::exp(std::complex<double>(0, m_wavenumber * minusPhi));
  }

  State m_state;
  double m_wavenumber;
  const double* m_maskData;
  const double* m_zernikesData;
  std::vector<double> m_alphas;
  Fourier::ComplexDft m_pupilToPsf;
  Fourier::ComplexDftBuffer& m_pupilAmplitude;
  Fourier::ComplexDftBuffer& m_psfAmplitude;
  Fourier::RealDftBuffer m_psfIntensity;
};

} // namespace Duffieux
} // namespace Phi

#endif // _PHIDUFFIEUX_MONOCHROMATICOPTICS_H
