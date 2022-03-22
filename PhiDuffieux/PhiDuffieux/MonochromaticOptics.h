// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef _PHIDUFFIEUX_MONOCHROMATICOPTICS_H
#define _PHIDUFFIEUX_MONOCHROMATICOPTICS_H

#include "PhiBox/StepperAlgo.h"
#include "PhiFourier/Dft.h"

#include <vector>

namespace Phi {
namespace Duffieux {

struct PupilAmplitude {
  using Prerequisite = void;
  using Return = const Fourier::ComplexDftBuffer&;
};

struct PsfAmplitude {
  using Prerequisite = PupilAmplitude;
  using Return = const Fourier::ComplexDftBuffer&;
};

struct PsfIntensity {
  using Prerequisite = PsfAmplitude;
  using Return = const Fourier::RealDftBuffer&;
};

/**
 * @brief Monochromatic data buffers and transforms with lazy evaluation.
 */
class MonochromaticOptics : public Framework::StepperAlgo<MonochromaticOptics> {
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
      m_wavenumber(2 * m_pi / lambda), m_maskData(mask.data()), m_zernikesData(basis.data()),
      m_alphas(std::move(alphas)), m_pupilToPsf(mask.shape()), m_pupilAmplitude(m_pupilToPsf.in()),
      m_psfAmplitude(m_pupilToPsf.out()), m_psfIntensity(mask.shape()) {}

  /**
   * @brief Update the wavelength.
   */
  void updateLambda(double lambda) {
    m_wavenumber = 2 * m_pi / lambda;
    reset();
  }

protected:
  /**
   * @brief Get the result of step `S`.
   */
  template <typename S>
  typename S::Return doGet() const;

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
    for (auto aIt = m_alphas.begin(); aIt != m_alphas.end(); ++aIt, ++zIt) {
      minusPhi -= (*aIt) * (*zIt);
    }
    return mask * std::exp(std::complex<double>(0, m_wavenumber * minusPhi));
  }

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
