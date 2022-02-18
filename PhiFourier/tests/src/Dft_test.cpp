// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "PhiFourier/Dft.h"

#include <boost/test/unit_test.hpp>

using namespace Phi::Fourier;

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Dft_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(dft_r2c_c2c_ic2c_ir2c_test) {

  constexpr long width = 5;
  constexpr long halfWidth = (width + 1) / 2;
  constexpr long height = 6;
  const Position shape {width, height};
  const Position halfShape {halfWidth, height};

  // Make r2c
  RealDft real({width, height});
  BOOST_TEST(real.in().shape() == shape);
  BOOST_TEST(real.in().owns());
  BOOST_TEST(real.in().data() != nullptr);
  BOOST_TEST(real.out().shape() == halfShape);
  BOOST_TEST(real.out().owns());
  BOOST_TEST(real.out().data() != nullptr);

  // Make ir2c
  auto inverseReal = real.inverse();
  BOOST_TEST(inverseReal.in().shape() == halfShape);
  BOOST_TEST(not inverseReal.in().owns());
  BOOST_TEST(inverseReal.in().data() == real.out().data());
  BOOST_TEST(inverseReal.out().shape() == shape);
  BOOST_TEST(not inverseReal.out().owns());
  BOOST_TEST(inverseReal.out().data() == real.in().data());

  // Make c2c
  auto complex = real.template compose<ComplexDft>(real.outShape());
  BOOST_TEST(complex.in().shape() == halfShape);
  BOOST_TEST(not complex.in().owns());
  BOOST_TEST(complex.in().data() == real.out().data());
  BOOST_TEST(complex.out().shape() == halfShape);
  BOOST_TEST(complex.out().owns());
  BOOST_TEST(complex.out().data() != nullptr);

  // Make ic2c
  auto inverseComplex = complex.inverse();
  BOOST_TEST(inverseComplex.in().shape() == halfShape);
  BOOST_TEST(not inverseComplex.in().owns());
  BOOST_TEST(inverseComplex.in().data() == complex.out().data());
  BOOST_TEST(inverseComplex.out().shape() == halfShape);
  BOOST_TEST(not inverseComplex.out().owns());
  BOOST_TEST(inverseComplex.out().data() == complex.in().data());

  // Fill r2c signal
  auto& signal = real.in();
  for (const auto& p : signal.domain()) {
    signal[p] = 1 + p[0] + p[1];
  }

  // Apply r2c and then c2c
  BOOST_TEST(real.in().data() != nullptr);
  BOOST_TEST(real.out().data() != nullptr);
  real.transform();
  BOOST_TEST(complex.in().data() != nullptr);
  BOOST_TEST(complex.out().data() != nullptr);
  complex.transform();

  // Apply ic2c and then ir2c
  BOOST_TEST(inverseComplex.in().data() != nullptr);
  BOOST_TEST(inverseComplex.out().data() != nullptr);
  inverseComplex.transform().normalize();
  BOOST_TEST(inverseReal.in().data() != nullptr);
  BOOST_TEST(inverseReal.out().data() != nullptr);
  inverseReal.transform().normalize();

  // Check values are recovered
  for (const auto& p : signal.domain()) {
    const auto expected = 1 + p[0] + p[1];
    const auto value = signal[p];
    BOOST_TEST(value > 0.99 * expected);
    BOOST_TEST(value < 1.01 * expected);
  }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
