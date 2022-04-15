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
  const std::vector<double> expected {4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 16};
  Sampling1D<const int> in(inData.size(), inData.data());
  Kernel1D<float> kernel(kernelData, 1);
  Sampling1D<double> out(outData.size(), outData.data());
  convolve1DZero(in, kernel, out);
  for (std::size_t i = 0; i < expected.size(); ++i) {
    BOOST_TEST(outData[i] == expected[i]);
  }
}

BOOST_AUTO_TEST_CASE(stepped_convolve1d_test) {
  const std::vector<int> inData {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
  const std::vector<float> kernelData {0.5, 1., 1.5, 2., 1.};
  std::vector<double> outData {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  const std::vector<double> expected {1, 2, 8.5, 4, 26, 6, 44, 8, 62, 10, 65, 12};
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

template <typename T>
std::vector<T> transposePad(const std::vector<T>& in, long stride) {
  std::vector<T> out(in.size() * stride, T());
  for (std::size_t i = 0; i < in.size(); ++i) {
    out[i * stride] = in[i];
  }
  return out;
}

BOOST_AUTO_TEST_CASE(transpose_pad_test) {
  const std::vector<int> in {1, 2, 3, 4};
  const auto out = transposePad(in, 3);
  const std::vector<int> expected {1, 0, 0, 2, 0, 0, 3, 0, 0, 4, 0, 0};
  BOOST_TEST(out == expected);
}

BOOST_AUTO_TEST_CASE(strided_convolve1d_test) {
  const auto inData = transposePad<int>({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14}, 2);
  const std::vector<float> kernelData {0.5, 1., 1.5, 2., 1.};
  auto outData = transposePad<double>({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}, 3);
  const auto expected = transposePad<double>({1, 2, 8.5, 4, 26, 6, 44, 8, 62, 10, 65, 12}, 3);
  Sampling1D<const int> in(inData.size() / 2, inData.data());
  in.from(1).step(3).stride(2);
  Kernel1D<float> kernel(kernelData, 3);
  Sampling1D<double> out(outData.size() / 3, outData.data());
  out.from(2).step(2).stride(3);
  convolve1DZero(in, kernel, out);
  for (std::size_t i = 0; i < expected.size(); ++i) {
    BOOST_TEST(outData[i] == expected[i]);
  }
}

BOOST_AUTO_TEST_CASE(standard_convolvexy_test) {
  VecRaster<int, 2> in({4, 3});
  std::fill(in.begin(), in.end(), 2);
  std::vector<int> kernelData {1, 1, 1};
  Kernel1D<int> kernel(kernelData, 1);
  VecRaster<int, 2> out(in.shape());
  convolveXyZero(in, kernel, out);
  const std::vector<int> expected {8, 12, 12, 8, 12, 18, 18, 12, 8, 12, 12, 8};
  BOOST_TEST(out.size() == expected.size());
  for (std::size_t i = 0; i < expected.size(); ++i) {
    BOOST_TEST(out.data()[i] == expected[i]);
  }
}

BOOST_AUTO_TEST_CASE(stepped_no_edge_convolvexy_test) {
  printf("\nSTEPPED NO EDGE\n\n");
  VecRaster<int, 2> in({4 * 3 + 2, 3 * 2 + 2});
  std::fill(in.begin(), in.end(), 2);
  std::vector<int> kernelData {1, 1, 1};
  Kernel1D<int> kernel(kernelData, 1);
  VecRaster<int, 2> out({4, 3});
  convolveXyZero(in, {{1, 1}, {4 * 3, 3 * 2}}, {3, 2}, kernel, out);
  const std::vector<int> expected(12, 18);
  BOOST_TEST(out.size() == expected.size());
  for (std::size_t i = 0; i < expected.size(); ++i) {
    BOOST_TEST(out.data()[i] == expected[i]);
  }
}

BOOST_AUTO_TEST_CASE(stepped_front_edge_convolvexy_test) {
  printf("\nSTEPPED FRONT EDGE\n\n");
  VecRaster<int, 2> in({4 * 3 + 2, 3 * 2 + 2});
  std::fill(in.begin(), in.end(), 2);
  std::vector<int> kernelData {1, 1, 1};
  Kernel1D<int> kernel(kernelData, 1);
  VecRaster<int, 2> out({4, 3});
  convolveXyZero(in, {{0, 0}, {3 * 3, 2 * 2}}, {3, 2}, kernel, out);
  const std::vector<int> expected {8, 12, 12, 12, 12, 18, 18, 18, 12, 18, 18, 18};
  BOOST_TEST(out.size() == expected.size());
  for (std::size_t i = 0; i < expected.size(); ++i) {
    BOOST_TEST(out.data()[i] == expected[i]);
  }
}

BOOST_AUTO_TEST_CASE(stepped_back_edge_convolvexy_test) {
  printf("\nSTEPPED BACK EDGE\n\n");
  VecRaster<int, 2> in({4 * 3 + 2, 3 * 2 + 2});
  std::fill(in.begin(), in.end(), 2);
  std::vector<int> kernelData {1, 1, 1};
  Kernel1D<int> kernel(kernelData, 1);
  VecRaster<int, 2> out({4, 3});
  convolveXyZero(in, {{4, 3}, {4 * 3 + 1, 3 * 2 + 1}}, {3, 2}, kernel, out);
  const std::vector<int> expected {18, 18, 18, 12, 18, 18, 18, 12, 12, 12, 12, 8};
  BOOST_TEST(out.size() == expected.size());
  for (std::size_t i = 0; i < expected.size(); ++i) {
    BOOST_TEST(out.data()[i] == expected[i]);
  }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
