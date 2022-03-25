// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef _PHIDUFFIEUX_BROADBANDSYSTEM_H
#define _PHIDUFFIEUX_BROADBANDSYSTEM_H

#include "PhiDuffieux/MonochromaticSystem.h"

namespace Phi {
namespace Duffieux {

class BroadbandSystem {

public:
  BroadbandSystem(std::vector<double> lambdas, MonochromaticOptics::Params optics, MonochromaticSystem system);

  template <typename S>
  typename S::Return get(long index) {
    return m_systems[index].template get<S>();
  }

private:
  std::vector<MonochromaticSystem> m_systems;
};

} // namespace Duffieux
} // namespace Phi

#endif // _PHIDUFFIEUX_BROADBANDSYSTEM_H
