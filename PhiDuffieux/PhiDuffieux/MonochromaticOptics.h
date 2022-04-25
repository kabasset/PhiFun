// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef _PHIDUFFIEUX_MONOCHROMATICOPTICS_H
#define _PHIDUFFIEUX_MONOCHROMATICOPTICS_H

#include "PhiBox/StepperPipeline.h"
#include "PhiFourier/Dft.h"

#include <algorithm>
#include <vector>

namespace Phi {
namespace Duffieux {

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
  struct Parameters {
    double wavelength; ///< The current wavelength
    Fourier::Position shape; ///< The input data shape
    const double* pupilMask; ///< The pupil mask data
    const double* zernikeBasis; ///< The Zernike basis data
    std::vector<double> zernikeCoefficients; ///< The Zernike coefficients
  };

  /**
   * @brief Constructor.
   */
  explicit MonochromaticOptics(Parameters parameters) :
      m_parameters(std::move(parameters)), m_wavenumber(2 * m_pi / m_parameters.wavelength),
      m_pupilToPsf(m_parameters.shape), m_psfIntensity(m_parameters.shape) {}

  /**
   * @brief Get the optical parameters.
   */
  const Parameters& parameters() const {
    return m_parameters;
  }

  /**
   * @brief Get the wavelength.
   */
  double wavelength() const {
    return m_parameters.wavelength;
  }

  /**
   * @brief Update the wavelength and optionally the Zernike coefficients.
   * @warning
   * This method resets the pipeline.
   */
  void update(double wavelength, const double* coefficients = nullptr) {
    m_parameters.wavelength = wavelength;
    m_wavenumber = 2 * m_pi / wavelength;
    if (coefficients) {
      std::copy_n(coefficients, m_parameters.zernikeCoefficients.size(), m_parameters.zernikeCoefficients.data());
    }
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
  inline std::complex<double> evalPhase(double mask, const double* zernikes) {
    double minusPhi = 0;
    auto zIt = zernikes;
    for (auto aIt = m_parameters.zernikeCoefficients.begin(); aIt != m_parameters.zernikeCoefficients.end();
         ++aIt, ++zIt) {
      minusPhi -= (*aIt) * (*zIt);
    } // Faster than inner_product
    const auto arg = m_wavenumber * minusPhi;
    return std::polar(mask, arg);
    // exp is slower when the argument is pure imaginary
    // See notes: https://en.cppreference.com/w/cpp/numeric/complex/exp
  }

  Parameters m_parameters; ///< The optical parameters
  double m_wavenumber; ///< The wave number
  Fourier::ComplexDft m_pupilToPsf; ///< The pupil amplitude to PSF amplitude DFT
  Fourier::RealDftBuffer m_psfIntensity; ///< The PSF intensity
};

/**
 * @brief Complex pupil amplitude computed from the phase and pupil mask.
 * @details
 * The phase is obtained from a Zernike basis and Zernike coefficients.
 */
struct PupilAmplitude : Framework::PipelineStep<void, const Fourier::ComplexDftBuffer&> {};

template <>
inline const Fourier::ComplexDftBuffer& MonochromaticOptics::doGet<PupilAmplitude>() {
  return m_pupilToPsf.in();
}

template <>
inline void MonochromaticOptics::doEvaluate<PupilAmplitude>() {
  const auto size = m_parameters.zernikeCoefficients.size();
  const double* maskIt = m_parameters.pupilMask;
  const double* zernikesIt = m_parameters.zernikeBasis;
  auto& amp = m_pupilToPsf.in();
  for (auto it = amp.begin(); it != amp.end(); ++it, ++maskIt, zernikesIt += size) {
    if (*maskIt != 0) {
      *it = evalPhase(*maskIt, zernikesIt);
    } else {
      *it = 0;
    }
  }
}

/**
 * @brief Complex PSF amplitude computed as the inverse DFT of the complex pupil amplitude.
 */
struct PsfAmplitude : Framework::PipelineStep<PupilAmplitude, const Fourier::ComplexDftBuffer&> {};

template <>
inline const Fourier::ComplexDftBuffer& MonochromaticOptics::doGet<PsfAmplitude>() {
  return m_pupilToPsf.out();
}

template <>
inline void MonochromaticOptics::doEvaluate<PsfAmplitude>() {
  m_pupilToPsf.transform().normalize();
}

/**
 * @brief PSF intensity computed as the norm squared of the PSF amplitude.
 */
struct PsfIntensity : Framework::PipelineStep<PsfAmplitude, const Fourier::RealDftBuffer&> {};

template <>
inline const Fourier::RealDftBuffer& MonochromaticOptics::doGet<PsfIntensity>() {
  return m_psfIntensity;
}

template <>
inline void MonochromaticOptics::doEvaluate<PsfIntensity>() {
  norm2(m_pupilToPsf.out(), m_psfIntensity);
}

} // namespace Duffieux
} // namespace Phi

#endif // _PHIDUFFIEUX_MONOCHROMATICOPTICS_H
