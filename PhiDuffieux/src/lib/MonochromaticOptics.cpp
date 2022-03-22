// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "PhiDuffieux/MonochromaticOptics.h"

namespace Phi {
namespace Duffieux {

template <>
const Fourier::ComplexDftBuffer& MonochromaticOptics::doGet<PupilAmplitude>() const {
  return m_pupilToPsf.in();
}

template <>
const Fourier::ComplexDftBuffer& MonochromaticOptics::doGet<PsfAmplitude>() const {
  return m_pupilToPsf.out();
}

template <>
const Fourier::RealDftBuffer& MonochromaticOptics::doGet<PsfIntensity>() const {
  return m_psfIntensity;
}

template <>
void MonochromaticOptics::doEvaluate<PupilAmplitude>() {
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
void MonochromaticOptics::doEvaluate<PsfAmplitude>() {
  m_pupilToPsf.transform().normalize();
}

template <>
void MonochromaticOptics::doEvaluate<PsfIntensity>() {
  norm2(m_pupilToPsf.out(), m_psfIntensity);
}

} // namespace Duffieux
} // namespace Phi
