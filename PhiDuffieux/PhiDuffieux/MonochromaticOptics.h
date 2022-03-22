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
   * @brief Optical parameters.
   */
  struct Params {
    double wavenumber;
    Fourier::Position shape;
    const double* maskData;
    const double* zernikesData;
    std::vector<double> alphas;
  };

  /**
   * @brief Constructor.
   * @param lambda The wavelength
   * @param diameter The pupil mask
   * @param basis The Zernike basis
   * @param alphas The Zernike coefficients
   */
  template <typename TMask, typename TZernikes>
  MonochromaticOptics(double lambda, const TMask& mask, const TZernikes& basis, std::vector<double> alphas) :
      m_params {2 * m_pi / lambda, mask.shape(), mask.data(), basis.data(), std::move(alphas)},
      m_pupilToPsf(m_params.shape), m_psfIntensity(m_params.shape) {}

  /**
   * @brief Update the wavelength.
   */
  void updateLambda(double lambda) {
    m_params.wavenumber = 2 * m_pi / lambda;
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
    for (auto aIt = m_params.alphas.begin(); aIt != m_params.alphas.end(); ++aIt, ++zIt) {
      minusPhi -= (*aIt) * (*zIt);
    }
    return mask * std::exp(std::complex<double>(0, m_params.wavenumber * minusPhi));
  }

  Params m_params;
  Fourier::ComplexDft m_pupilToPsf;
  Fourier::RealDftBuffer m_psfIntensity;
};

} // namespace Duffieux
} // namespace Phi

#endif // _PHIDUFFIEUX_MONOCHROMATICOPTICS_H
