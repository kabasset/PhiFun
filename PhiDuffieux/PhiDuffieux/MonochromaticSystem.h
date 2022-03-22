// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef _PHIDUFFIEUX_MONOCHROMATICSYSTEM_H
#define _PHIDUFFIEUX_MONOCHROMATICSYSTEM_H

#include "PhiBox/ImageProcessing.h"
#include "PhiDuffieux/MonochromaticOptics.h"

namespace Phi {
namespace Duffieux {

class MonochromaticSystem {

public:
  MonochromaticSystem(MonochromaticOptics& optics, long side) :
      m_optics(optics), m_psf(m_optics.m_psfIntensity), m_psfToTf(m_psf.shape(), m_psf.data()), m_tfToPsf({side, side}),
      m_stf(m_psfToTf.out()) {}

  const Fourier::ComplexDftBuffer& evalOpticalTf() {
    m_psfToTf.transform();
    return m_stf;
  }

  const Fourier::ComplexDftBuffer& evalSystemTf(const Fourier::ComplexDftBuffer& nonOpticalTf) {
    m_stf *= nonOpticalTf;
    return m_stf;
  }

  /**
   * @brief Warp the system transfer function.
   * @details
   * The distortion model is linear as follows:
   * - `u = ux * x + uy * y`,
   * - `v = vx * x + vy * y`;
   * 
   * where the coordinates are normalized and coefficients depend linearly on the wavelength.
   */
  const Fourier::ComplexDftBuffer& warpSystemTf(double ux, double uy, double vx, double vy) {
    const auto width = m_tfToPsf.in().shape()[0];
    const auto height = m_tfToPsf.in().shape()[1];
    const double xFactor = double(m_stf.shape()[0] - 1) / (width - 1);
    const double yFactor = double(m_stf.shape()[1] - 1) / (height - 1);
    const double uxFactor = ux * xFactor;
    const double uyFactor = uy * yFactor;
    const double vxFactor = vx * xFactor;
    const double vyFactor = vy * yFactor;
    auto* it = m_tfToPsf.in().begin();
    for (long y = 0; y < height; ++y) {
      const double uyValue = uyFactor * y;
      const double vyValue = vyFactor * y;
      for (long x = 0; x < width; ++x, ++it) {
        const double u = uxFactor * x + uyValue;
        const double v = vxFactor * x + vyValue;
        *it = Image2D::bilinear<std::complex<double>>(m_stf, u, v);
      }
    }
    return m_tfToPsf.in();
  }

  const Fourier::RealDftBuffer& evalSystemPsf() {
    m_tfToPsf.transform().normalize();
    return m_tfToPsf.out();
  }

private:
  MonochromaticOptics& m_optics;
  Fourier::RealDftBuffer& m_psf;
  Fourier::RealDft m_psfToTf;
  Fourier::RealDft::Inverse m_tfToPsf;
  Fourier::ComplexDftBuffer& m_stf;
};

} // namespace Duffieux
} // namespace Phi

#endif // _PHIDUFFIEUX_MONOCHROMATICSYSTEM_H
