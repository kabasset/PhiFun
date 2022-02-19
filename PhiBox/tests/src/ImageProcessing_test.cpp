// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "EleFits/SifFile.h" // FIXME rm
#include "PhiBox/ImageProcessing.h"

#include <boost/test/unit_test.hpp>

using namespace Phi::Image;

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(ImageProcessing_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(bilinear_test) {
  const Euclid::Fits::VecRaster<int> input({2, 2}, std::vector<int> {0, 2, 3, 1});
  Euclid::Fits::VecRaster<double> output({100, 50});
  const double xFactor = 1. / (output.length<0>() - 1);
  const double yFactor = 1. / (output.length<1>() - 1);
  for (const auto& p : output.domain()) {
    output[p] = bilinear<double>(input, p[0] * xFactor, p[1] * yFactor);
  }

  // FIXME Strict equality?
  BOOST_TEST((output[{0, 0}]) == 0);
  BOOST_TEST((output.at({-1, 0})) == 2);
  BOOST_TEST((output.at({0, -1})) == 3);
  BOOST_TEST((output.at({-1, -1})) == 1);

  Euclid::Fits::SifFile f("/tmp/bilinear.fits", Euclid::Fits::FileMode::Overwrite);
  f.writeRaster(output);
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
