// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "PhiBox/StepperPipeline.h"

#include <boost/test/unit_test.hpp>

//-----------------------------------------------------------------------------

using namespace Phi::Framework;

struct Step0 : PipelineStep<void, char> {
  Return value = 0;
};

struct Step1a : PipelineStep<Step0, short> {
  Return value = 0;
};

struct Step1b : PipelineStep<Step0, int> {
  Return value = 0;
};

struct Step2 : PipelineStep<std::tuple<Step1a, Step1b>, long> {
  Return value = 0;
};

class Dag : public StepperPipeline<Dag> {
public:
  Step0::Return get0() const {
    return m_0.value;
  }

  Step1a::Return get1a() const {
    return m_1a.value;
  }

  Step1b::Return get1b() const {
    return m_1b.value;
  }

  Step2::Return get2() const {
    return m_2.value;
  }

protected:
  template <typename S>
  void doEvaluate();

  template <typename S>
  typename S::Return doGet();

private:
  char m_value = 0;
  Step0 m_0;
  Step1a m_1a;
  Step1b m_1b;
  Step2 m_2;
};

template <>
void Dag::doEvaluate<Step0>() {
  ++m_value;
  m_0.value = m_value;
}

template <>
void Dag::doEvaluate<Step1a>() {
  ++m_value;
  m_1a.value = m_value;
}

template <>
void Dag::doEvaluate<Step1b>() {
  ++m_value;
  m_1b.value = m_value;
}

template <>
void Dag::doEvaluate<Step2>() {
  ++m_value;
  m_2.value = m_value;
}

template <>
Step0::Return Dag::doGet<Step0>() {
  return m_0.value;
}

template <>
Step1a::Return Dag::doGet<Step1a>() {
  return m_1a.value;
}

template <>
Step1b::Return Dag::doGet<Step1b>() {
  return m_1b.value;
}

template <>
Step2::Return Dag::doGet<Step2>() {
  return m_2.value;
}

BOOST_AUTO_TEST_SUITE(StepperPipeline_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(back_and_forth_test) {

  // Init
  Dag dag;
  BOOST_TEST(dag.get0() == 0);
  BOOST_TEST(dag.get1a() == 0);
  BOOST_TEST(dag.get1b() == 0);
  BOOST_TEST(dag.get2() == 0);
  BOOST_TEST(dag.milliseconds<Step0>() < 0);
  BOOST_TEST(dag.milliseconds<Step1a>() < 0);
  BOOST_TEST(dag.milliseconds<Step1b>() < 0);
  BOOST_TEST(dag.milliseconds<Step2>() < 0);
  const auto a = dag.get<Step1a>();
  BOOST_TEST(dag.get0() == 1);
  BOOST_TEST(dag.get1a() == 2);
  BOOST_TEST(dag.get1b() == 0);
  BOOST_TEST(dag.get2() == 0);
  BOOST_TEST(a == 2);
  BOOST_TEST(dag.milliseconds<Step0>() > 0);
  BOOST_TEST(dag.milliseconds<Step1a>() > 0);
  BOOST_TEST(dag.milliseconds<Step1b>() < 0);
  BOOST_TEST(dag.milliseconds<Step2>() < 0);

  // Back
  const auto o = dag.get<Step0>();
  BOOST_TEST(dag.get0() == 1);
  BOOST_TEST(dag.get1a() == 2);
  BOOST_TEST(dag.get1b() == 0);
  BOOST_TEST(dag.get2() == 0);
  BOOST_TEST(o == 1);
  BOOST_TEST(dag.milliseconds<Step0>() > 0);
  BOOST_TEST(dag.milliseconds<Step1a>() > 0);
  BOOST_TEST(dag.milliseconds<Step1b>() < 0);
  BOOST_TEST(dag.milliseconds<Step2>() < 0);

  // Forth
  const auto z = dag.get<Step2>();
  BOOST_TEST(dag.get0() == 1);
  BOOST_TEST(dag.get1a() == 2);
  BOOST_TEST(dag.get1b() == 3);
  BOOST_TEST(dag.get2() == 4);
  BOOST_TEST(z == 4);
  BOOST_TEST(dag.milliseconds<Step0>() > 0);
  BOOST_TEST(dag.milliseconds<Step1a>() > 0);
  BOOST_TEST(dag.milliseconds<Step1b>() > 0);
  BOOST_TEST(dag.milliseconds<Step2>() > 0);
}

BOOST_AUTO_TEST_CASE(all_in_one_test) {
  Dag dag;
  const auto z = dag.get<Step2>();
  BOOST_TEST(dag.get0() == 1);
  BOOST_TEST(dag.get1a() == 2);
  BOOST_TEST(dag.get1b() == 3);
  BOOST_TEST(dag.get2() == 4);
  BOOST_TEST(z == 4);
  BOOST_TEST(dag.milliseconds<Step0>() > 0);
  BOOST_TEST(dag.milliseconds<Step1a>() > 0);
  BOOST_TEST(dag.milliseconds<Step1b>() > 0);
  BOOST_TEST(dag.milliseconds<Step2>() > 0);
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
