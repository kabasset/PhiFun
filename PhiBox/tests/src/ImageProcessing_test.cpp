// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "EleFits/SifFile.h" // FIXME rm
#include "PhiBox/ImageProcessing.h"

#include <boost/test/unit_test.hpp>

using namespace Phi::Image2D;

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(ImageProcessing_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(bilinear_test) {
  const Euclid::Fits::VecRaster<int> input({2, 2}, std::vector<int> {0, 2, 3, 1});
  Euclid::Fits::VecRaster<double> output({100, 50});
  const double xFactor = 1. / (output.length<0>() - 1);
  const double yFactor = 1. / (output.length<1>() - 1);
  for (const auto& p : output.domain()) {
    output[p] = bilinearNearest<double>(input, p[0] * xFactor, p[1] * yFactor);
  }

  // FIXME Approx instead of floor
  BOOST_TEST((int(output[{0, 0}] + .5)) == 0);
  BOOST_TEST((int(output.at({-1, 0}) + .5)) == 2);
  BOOST_TEST((int(output.at({0, -1}) + .5)) == 3);
  BOOST_TEST((int(output.at({-1, -1}) + .5)) == 1);

  Euclid::Fits::SifFile f("/tmp/bilinear.fits", Euclid::Fits::FileMode::Overwrite); // FIXME rm
  f.writeRaster(output);
}

BOOST_AUTO_TEST_CASE(sampling1d_test) {
  const std::vector<long> inData {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  Sampling1D<const long> in(inData.size(), inData.data());
  in.from(1).to(6).step(2);
  const std::vector<long> expected {2, 4, 6};
  BOOST_TEST(in.count() == expected.size());
  BOOST_TEST(*in.begin() == 2);
  BOOST_TEST(*in.end() == 8);
  long i = 0;
  for (auto v : in) {
    BOOST_TEST(v == expected[i]);
    ++i;
  }
  BOOST_TEST(i == expected.size());
}

BOOST_AUTO_TEST_CASE(standard_convolve1d_test) {
  const std::vector<int> inData {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
  const std::vector<float> kernelData {0.5, 1., 1.5};
  std::vector<double> outData(inData.size(), 0.);
  std::vector<double> expected {4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 16};
  Sampling1D<const int> in(inData.size(), inData.data());
  Kernel1D<float> kernel(kernelData, 1);
  Sampling1D<double> out(outData.size(), outData.data());
  convolve1DZero(in, kernel, out);
  for (std::size_t i = 0; i < expected.size(); ++i) {
    BOOST_TEST(outData[i] == expected[i]);
  }
}

BOOST_AUTO_TEST_CASE(steped_convolve1d_test) {
  const std::vector<int> inData {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
  const std::vector<float> kernelData {0.5, 1., 1.5, 2., 1.};
  std::vector<double> outData {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  std::vector<double> expected {1, 2, 8.5, 4, 26, 6, 44, 8, 62, 10, 65, 12};
  Sampling1D<const int> in(inData.size(), inData.data());
  in.from(1);
  in.step(3);
  Kernel1D<float> kernel(kernelData, 3);
  Sampling1D<double> out(outData.size(), outData.data());
  out.from(2);
  out.step(2);
  convolve1DZero(in, kernel, out);
  for (std::size_t i = 0; i < expected.size(); ++i) {
    BOOST_TEST(outData[i] == expected[i]);
  }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
