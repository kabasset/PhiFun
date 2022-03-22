// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef _PHIDUFFIEUX_BROADBANDSYSTEM_H
#define _PHIDUFFIEUX_BROADBANDSYSTEM_H

#include "PhiDuffieux/MonochromaticSystem.h"

#include <vector>

namespace Phi {
namespace Duffieux {

class BroadbandSystem {

public:
  BroadbandSystem(std::vector<double> lambdas);

  MonochromaticSystem& system(long index) {
    return m_systems[index];
  }

  MonochromaticOptics& optics(long index) {
    return m_optics[index];
  }

private:
  std::vector<MonochromaticOptics> m_optics;
  std::vector<MonochromaticSystem> m_systems;
};

} // namespace Duffieux
} // namespace Phi

#endif // _PHIDUFFIEUX_BROADBANDSYSTEM_H
