// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "PhiBox/StepperAlgo.h"

#include <boost/test/unit_test.hpp>

//-----------------------------------------------------------------------------

struct StepA {
  using Prerequisite = void;
  using Return = int;
  Return value = 0;
};

struct StepB {
  using Prerequisite = StepA;
  using Return = char;
  Return value = 'z';
};

class ABAlgo : public Phi::Framework::StepperAlgo<ABAlgo> {
public:
  StepA::Return getA() const {
    return m_a.value;
  }

  StepB::Return getB() const {
    return m_b.value;
  }

protected:
  template <typename S>
  void doEvaluate();

  template <typename S>
  typename S::Return doGet();

private:
  StepA m_a;
  StepB m_b;
};

template <>
void ABAlgo::doEvaluate<StepA>() {
  ++m_a.value;
}

template <>
void ABAlgo::doEvaluate<StepB>() {
  --m_b.value;
}

template <>
StepA::Return ABAlgo::doGet<StepA>() {
  return m_a.value;
}

template <>
StepB::Return ABAlgo::doGet<StepB>() {
  return m_b.value;
}

BOOST_AUTO_TEST_SUITE(StepperAlgo_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(caching_test) {
  ABAlgo algo;
  BOOST_TEST(algo.getA() == 0);
  BOOST_TEST(algo.getB() == 'z');
  const auto b = algo.get<StepB>();
  BOOST_TEST(algo.getA() == 1);
  BOOST_TEST(algo.getB() == 'y');
  BOOST_TEST(b == 'y');
  const auto a = algo.get<StepA>();
  BOOST_TEST(algo.getA() == 1);
  BOOST_TEST(algo.getB() == 'y');
  BOOST_TEST(a == 1);
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
