// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "PhiFourier/DftBuffer.h"

#include <boost/test/unit_test.hpp>

using namespace Phi::Fourier;

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(DftBuffer_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(owned_and_shared_test) {
  constexpr long width = 3;
  constexpr long height = 4;
  constexpr long size = width * height;
  DftBuffer<double> owner({width, height});
  BOOST_TEST(owner.size() == size);
  const DftBuffer<double> sharer({width, height}, owner.data());
  BOOST_TEST(not sharer.owns());
  BOOST_TEST(sharer.size() == size);
  BOOST_TEST(sharer.data() == owner.data());
  DftBuffer<const double> observer({width, height}, owner.data());
  BOOST_TEST(observer.size() == size);
  BOOST_TEST(observer.data() == owner.data());
  for (auto& v : owner) {
    v = 0;
  }
  for (const auto& p : owner.domain()) {
    BOOST_TEST(owner[p] == 0);
    BOOST_TEST(sharer[p] == 0);
    BOOST_TEST(observer[p] == 0);
  }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
