// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef _PHIDUFFIEUX_MONOCHROMATICSYSTEM_H
#define _PHIDUFFIEUX_MONOCHROMATICSYSTEM_H

#include "PhiDuffieux/MonochromaticOptics.h"

namespace Phi {
namespace Duffieux {

class MonochromaticSystem {

public:
  MonochromaticSystem(MonochromaticOptics& optics) :
      m_optics(optics), m_psf(m_optics.m_psfIntensity), m_psfToTf(m_psf.shape(), m_psf.data()),
      m_tfToPsf(m_psfToTf.inverse()), m_stf(m_psfToTf.out()) {}

  const Fourier::ComplexDftBuffer& evalOpticalTf() {
    m_psfToTf.transform();
    return m_stf;
  }

  const Fourier::ComplexDftBuffer& evalSystemTf(const Fourier::ComplexDftBuffer& nonOpticalTf) {
    m_stf.apply(
        [](const auto optical, const auto others) {
          return optical * others;
        },
        nonOpticalTf);
    return m_stf;
  }

  const Fourier::RealDftBuffer& evalSystemPsf() {
    m_tfToPsf.transform();
    return m_psf;
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
