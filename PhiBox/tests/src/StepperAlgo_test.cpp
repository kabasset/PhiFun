// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "PhiBox/StepperAlgo.h"

#include <boost/test/unit_test.hpp>

//-----------------------------------------------------------------------------

struct StepA {
  using Prerequisite = void;
  using Value = int;
  Value value = 0;
};

struct StepB {
  using Prerequisite = StepA;
  using Value = char;
  Value value = 'z';
};

class ABAlgo : public Phi::Framework::StepperAlgo<ABAlgo> {
public:
  StepA::Value getA() const {
    return m_a.value;
  }

  StepB::Value getB() const {
    return m_b.value;
  }

protected:
  template <typename S>
  void doEvaluate();

  template <typename S>
  typename S::Value doGet() const;

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
StepA::Value ABAlgo::doGet<StepA>() const {
  return m_a.value;
}

template <>
StepB::Value ABAlgo::doGet<StepB>() const {
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
