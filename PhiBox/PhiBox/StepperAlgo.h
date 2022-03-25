// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef _PHIBOX_STEPPERALGO_H
#define _PHIBOX_STEPPERALGO_H

#include <set>
#include <tuple>
#include <typeindex>
#include <utility>

namespace Phi {
namespace Framework {

/**
 * @brief Traits class which gives the cardinality (number of elemonts) of a type.
 * @details
 * Cardinality of:
 * - `void` is 0;
 * - a tuple is its size;
 * - any other type is 1.
 */
template <typename T>
struct TypeCardinality {
  static constexpr std::size_t value = 1;
};

template <>
struct TypeCardinality<void> {
  static constexpr std::size_t value = 0;
};

template <typename... Ts>
struct TypeCardinality<std::tuple<Ts...>> {
  static constexpr std::size_t value = sizeof...(Ts);
};

template <typename S>
constexpr int prerequisiteCardinality() {
  return TypeCardinality<typename S::Prerequisite>::value;
}

/**
 * @brief A pipeline or directed acyclic graph (DAG) which can be run step-by-step using lazy evaluation.
 * @details
 * The only public method, `get<S>()`, returns the value of step `S`.
 * If not already done, the prerequisites of `S` are first triggered.
 * 
 * This class relies on the CRTP.
 * Child classes should provide a specialization of the following methods for each step `S`:
 * - `void doEvaluate<S>()`, which evaluates `S` assuming upstream tasks were already computed;
 * - `S::Return doGet<S>()`, which returns the computed value of `S`.
 * 
 * A step `S` is a class which contains the following type definitions:
 * - `Return` is the return value type of `get<S>`;
 * - `Prerequisite` is (are) the step(s) which must be run prior to `S`, or `void` if there is no prerequisite;
 *   Multiple prerequisites are describled with tuples.
 */
template <typename Derived>
class StepperAlgo {

public:
  /**
   * @brief Evaluation of step `S` which has no prerequisite.
   */
  template <typename S>
  std::enable_if_t<prerequisiteCardinality<S>() == 0, typename S::Return> get() {
    return evaluateGet<S>();
  }

  /**
   * @brief Lazy evaluation of step `S` which triggers a single prerequisite evaluation if needed.
   */
  template <typename S>
  std::enable_if_t<prerequisiteCardinality<S>() == 1, typename S::Return> get() {
    get<typename S::Prerequisite>();
    return evaluateGet<S>();
  }

  /**
   * @brief Lazy evaluation of step `S` which triggers multiple prerequisite evaluations if needed.
   */
  template <typename S>
  std::enable_if_t<prerequisiteCardinality<S>() >= 2, typename S::Return> get() {
    getMultiple<typename S::Prerequisite>(std::make_index_sequence<prerequisiteCardinality<S>()> {});
    return evaluateGet<S>();
  }

protected:
  /**
   * @brief Reset to initial step.
   */
  void reset() {
    m_done.clear();
  }

private:
  template <typename STuple, std::size_t... Is>
  void getMultiple(std::index_sequence<Is...>) {
    using mockUnpack = int[];
    (void)mockUnpack {0, (get<std::tuple_element_t<Is, STuple>>(), void(), 0)...};
  }

  /**
   * @brief Access to protected methods of `Derived`.
   */
  template <typename S>
  struct Accessor : Derived {

    /**
     * @brief Call `algo.doEvaluate<S>()`.
     */
    static void evaluate(Derived& algo) {
      auto f = &Accessor::template doEvaluate<S>;
      (algo.*f)();
    }

    /**
     * @brief Call `algo.doGet<S>()`.
     */
    static typename S::Return get(Derived& algo) {
      auto f = &Accessor::template doGet<S>;
      return (algo.*f)();
    }
  };

  /**
   * @brief Run step `S` if not done.
   */
  template <typename S>
  typename S::Return evaluateGet() {
    const auto index = std::type_index(typeid(S));
    if (m_done.count(index) == 0) { // FIXME not thread-safe
      Accessor<S>::evaluate(derived());
      m_done.emplace(index);
    }
    return Accessor<S>::get(derived());
  }

  /**
   * @brief Cast as `Derived`.
   */
  inline Derived& derived() {
    return static_cast<Derived&>(*this);
  }

private:
  /**
   * @brief The set of performed steps.
   */
  std::set<std::type_index> m_done;
};

} // namespace Framework
} // namespace Phi

#endif // _PHIBOX_STEPPERALGO_H
