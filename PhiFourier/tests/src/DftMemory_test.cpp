// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "PhiFourier/Dft.h"
#include "PhiFourier/DftBuffer.h"
#include "PhiFourier/DftMemory.h"

#include <boost/test/unit_test.hpp>

using namespace Phi::Fourier;

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(DftMemory_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(allocate_plan_test) {
  const Position shape {6, 8};
  const Position halfShape {3, 8};
  RealDftBuffer rin(shape);
  RealDftBuffer rout(shape);
  ComplexDftBuffer cin(shape);
  ComplexDftBuffer cout(shape);
  auto rc = FftwAllocator::createPlan<RealDftType>(rin, cout);
  BOOST_TEST(rc != nullptr);
  FftwAllocator::destroyPlan(rc);
  auto irc = FftwAllocator::createPlan<Inverse<RealDftType>>(cout, rin);
  BOOST_TEST(irc != nullptr);
  FftwAllocator::destroyPlan(irc);
  auto cc = FftwAllocator::createPlan<ComplexDftType>(cin, cout);
  BOOST_TEST(cc != nullptr);
  FftwAllocator::destroyPlan(cc);
  auto icc = FftwAllocator::createPlan<Inverse<ComplexDftType>>(cout, cin);
  BOOST_TEST(icc != nullptr);
  FftwAllocator::destroyPlan(icc);
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
