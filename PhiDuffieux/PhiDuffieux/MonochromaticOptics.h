// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef _PHIDUFFIEUX_MONOCHROMATICOPTICS_H
#define _PHIDUFFIEUX_MONOCHROMATICOPTICS_H

#include "PhiBox/StepperPipeline.h"
#include "PhiFourier/Dft.h"

#include <vector>

namespace Phi {
namespace Duffieux {

/**
 * @brief Complex pupil amplitude computed from the phase and pupil mask.
 * @details
 * The phase is obtained from a Zernike basis and Zernike coefficients.
 */
struct PupilAmplitude {
  using Prerequisite = void;
  using Return = const Fourier::ComplexDftBuffer&;
};

/**
 * @brief Complex PSF amplitude computed as the inverse DFT of the complex pupil amplitude.
 */
struct PsfAmplitude {
  using Prerequisite = PupilAmplitude;
  using Return = const Fourier::ComplexDftBuffer&;
};

/**
 * @brief PSF intensity computed as the norm squared of the PSF amplitude.
 */
struct PsfIntensity {
  using Prerequisite = PsfAmplitude;
  using Return = const Fourier::RealDftBuffer&;
};

/**
 * @brief Monochromatic optical model.
 * @details
 * Contains parameters and data buffers and provides step-by-step transforms.
 * @see StepperPipeline
 */
class MonochromaticOptics : public Framework::StepperPipeline<MonochromaticOptics> {
  friend class MonochromaticSystem;

  /**
   * @brief 3.14 ;)
   */
  static constexpr double m_pi = std::acos(-1.);

public:
  /**
   * @brief Optical parameters.
   */
  struct Params {

    /**
     * @brief Constructor.
     * @param lambda The wavelength
     * @param diameter The pupil mask
     * @param basis The Zernike basis
     * @param alphas The Zernike coefficients
     */
    template <typename TMask, typename TZernikes>
    Params(double lambda, const TMask& mask, const TZernikes& basis, std::vector<double> coefficients) :
        wavelength(lambda), wavenumber(2 * m_pi / lambda), shape(mask.shape()), maskData(mask.data()),
        zernikesData(basis.data()), alphas(std::move(coefficients)) {}

    double wavelength; ///< The current wavelength
    double wavenumber; ///< The current wave number
    Fourier::Position shape; ///< The logical data shape
    const double* maskData; ///< The pupil mask data
    const double* zernikesData; ///< The Zernike basis data
    std::vector<double> alphas; ///< The Zernike coefficients
  };

  /**
   * @brief Constructor.
   */
  explicit MonochromaticOptics(Params params) :
      m_params(std::move(params)), m_pupilToPsf(m_params.shape), m_psfIntensity(m_params.shape) {}

  /**
   * @brief Update the wavelength.
   */
  void updateLambda(double lambda) {
    m_params.wavelength = lambda;
    m_params.wavenumber = 2 * m_pi / lambda;
    reset();
  }

protected:
  /**
   * @brief Get the result of step `S`.
   */
  template <typename S>
  typename S::Return doGet();

  /**
   * @brief Evaluate step `S`.
   */
  template <typename S>
  void doEvaluate();

private:
  /**
   * @brief Evaluate some local phase.
   * @param mask The local pupil mask value
   * @param zernikes The local Zernike basis values
   */
  std::complex<double> evalPhase(double mask, const double* zernikes) {
    double minusPhi = 0;
    auto zIt = zernikes;
    for (auto aIt = m_params.alphas.begin(); aIt != m_params.alphas.end(); ++aIt, ++zIt) {
      minusPhi -= (*aIt) * (*zIt);
    }
    return mask * std::exp(std::complex<double>(0, m_params.wavenumber * minusPhi));
  }

  Params m_params; ///< The optical parameters
  Fourier::ComplexDft m_pupilToPsf; ///< The pupil amplitude to PSF amplitude DFT
  Fourier::RealDftBuffer m_psfIntensity; ///< The PSF intensity
};

template <>
inline const Fourier::ComplexDftBuffer& MonochromaticOptics::doGet<PupilAmplitude>() {
  return m_pupilToPsf.in();
}

template <>
inline const Fourier::ComplexDftBuffer& MonochromaticOptics::doGet<PsfAmplitude>() {
  return m_pupilToPsf.out();
}

template <>
inline const Fourier::RealDftBuffer& MonochromaticOptics::doGet<PsfIntensity>() {
  return m_psfIntensity;
}

template <>
inline void MonochromaticOptics::doEvaluate<PupilAmplitude>() {
  const auto size = m_params.alphas.size();
  const double* maskIt = m_params.maskData;
  const double* zernikesIt = m_params.zernikesData;
  auto& amp = m_pupilToPsf.in();
  for (auto it = amp.begin(); it != amp.end(); ++it, ++maskIt, zernikesIt += size) {
    if (*maskIt != 0) {
      *it = evalPhase(*maskIt, zernikesIt);
    } else {
      *it = 0;
    }
  }
}

template <>
inline void MonochromaticOptics::doEvaluate<PsfAmplitude>() {
  m_pupilToPsf.transform().normalize();
}

template <>
inline void MonochromaticOptics::doEvaluate<PsfIntensity>() {
  norm2(m_pupilToPsf.out(), m_psfIntensity);
}

} // namespace Duffieux
} // namespace Phi

#endif // _PHIDUFFIEUX_MONOCHROMATICOPTICS_H
