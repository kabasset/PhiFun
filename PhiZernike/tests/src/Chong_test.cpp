// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "EleFits/SifFile.h" // FIXME rm
#include "EleFitsValidation/Chronometer.h" // FIXME rm
#include "PhiZernike/Chong.h"

#include <boost/test/unit_test.hpp>

using namespace Phi::Zernike;

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Chong_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(save_zernikes_test) {
  Euclid::Fits::Validation::Chronometer<std::chrono::milliseconds> chrono;
  constexpr long side = 1024;
  constexpr long orders = 10;
  chrono.start();
  const auto zxy = chongBasis(side, orders);
  chrono.stop();
  printf("Elapsed: %lims\n", chrono.last().count());
  Euclid::Fits::VecRaster<double, 3> xyz({side, side, zxy.shape()[0]});
  for (const auto& p : xyz.domain()) {
    xyz[p] = zxy[{p[2], p[1], p[0]}];
  }
  Euclid::Fits::SifFile f("/tmp/zernike.fits", Euclid::Fits::FileMode::Overwrite);
  f.writeRaster(xyz);
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
