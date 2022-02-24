// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "PhiBox/SplineIntegrator.h"

#include <boost/test/unit_test.hpp>
#include <numeric> // inner_product

using namespace Phi::Spline;

void testApprox(double value, double ref, double atol = 0.0001) {
  BOOST_TEST(std::abs(value - ref) < atol);
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(SplineIntegrator_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(constant_interpolation_test) {
  const std::vector<double> u {0., 0.5, 1., 2., 4.};
  const std::vector<double> x {0.0, 0.2, 1.5, 2.0, 2.3, 2.7, 3.0, 3.9};
  SplineIntegrator integrator(u, x);
  const std::vector<double> y(u.size(), 3.14);
  const auto z = integrator.knotZ(y.data());
  const auto s = integrator.interpolate(y.data(), z.data());
  BOOST_TEST(s.size() == x.size());
  for (const auto v : s) {
    testApprox(v, 3.14);
  }
}

BOOST_AUTO_TEST_CASE(affine_interpolation_test) {
  const std::vector<double> u {0., 0.5, 1., 2., 4.};
  const std::vector<double> x {0.0, 0.2, 1.5, 2.0, 2.3, 2.7, 3.0, 3.9};
  SplineIntegrator integrator(u, x);
  std::vector<double> y(u.size());
  std::transform(u.begin(), u.end(), y.begin(), [](auto v) {
    return 2. * v + 1.;
  });
  const auto z = integrator.knotZ(y.data());
  const auto s = integrator.interpolate(y.data(), z.data());
  BOOST_TEST(s.size() == x.size());
  for (std::size_t i = 0; i < x.size(); ++i) {
    testApprox(s[i], 2. * x[i] + 1);
  }
}

BOOST_AUTO_TEST_CASE(sin_interpolation_test) {
  const std::vector<double> u {0, 1, 2.5, 3.6, 5, 7, 8.1, 10};
  std::vector<double> x(40);
  for (std::size_t i = 0; i < x.size(); ++i) {
    x[i] = .25 * i;
  }
  SplineIntegrator integrator(u, x);
  std::vector<double> y(u.size());
  std::transform(u.begin(), u.end(), y.begin(), [](auto v) {
    return std::sin(v);
  });
  const auto z = integrator.knotZ(y.data());
  const auto s = integrator.interpolate(y.data(), z.data());
  BOOST_TEST(s.size() == x.size());
  for (std::size_t i = 0; i < x.size(); ++i) {
    testApprox(s[i], std::sin(x[i]), .1);
  }
}

BOOST_AUTO_TEST_CASE(linear_integration_test) {
  const std::vector<double> u {0., 0.5, 1., 2., 4.};
  const std::vector<double> x {0.0, 0.2, 1.5, 2.0, 2.3, 2.7, 3.0, 3.9};
  SplineIntegrator integrator(u, x);
  const std::vector<double> y(u.size(), 3.14);
  const auto z = integrator.knotZ(y.data());
  BOOST_TEST(z.size() == u.size());

  const auto sum = integrator.integrate(y.data(), z.data(), x.data());
  const auto expected = (0.2 + 1.5 + 2.0 + 2.3 + 2.7 + 3.0 + 3.9) * 3.14;
  testApprox(sum, expected);

  const auto s = integrator.interpolate(y.data(), z.data());
  BOOST_TEST(s.size() == x.size());
  const auto fromS = std::inner_product(s.begin(), s.end(), x.begin(), 0.);
  testApprox(fromS, expected);
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
