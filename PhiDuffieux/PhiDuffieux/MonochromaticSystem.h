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
  using Prerequisite = void;
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
class MonochromaticSystem : public Framework::StepperAlgo<MonochromaticSystem> {

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
  MonochromaticOptics& optics() {
    return m_optics;
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
  Params m_params; ///< The system parameters
  MonochromaticOptics m_optics; ///< The optical model
  Fourier::RealDft m_psfToTf; ///< The unwarped PSF to TF transform
  Fourier::RealDft::Inverse m_tfToPsf; ///< The warped TF to PSF transform
};

} // namespace Duffieux
} // namespace Phi

#endif // _PHIDUFFIEUX_MONOCHROMATICSYSTEM_H
