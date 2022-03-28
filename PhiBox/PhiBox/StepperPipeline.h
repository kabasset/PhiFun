// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef _PHIBOX_STEPPERPIPELINE_H
#define _PHIBOX_STEPPERPIPELINE_H

#include <chrono>
#include <map>
#include <tuple>
#include <typeindex>
#include <utility>

namespace Phi {
namespace Framework {

/**
 * @brief Traits class which gives the cardinality (number of elements) of a type.
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

/**
 * @brief Cardinality of a step's prerequisite.
 * @see `TypeCardinality`
 */
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
class StepperPipeline {

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

  /**
   * @brief Check whether some step `S` has already been evaluated.
   */
  template <typename S>
  bool evaluated() const {
    return m_milliseconds.find(key<S>()) != m_milliseconds.end();
  }

  /**
   * @brief Get the elapsed time of step `S`.
   * @return The time in millisecond if the step was evaluated, or -1 otherwise.
   */
  template <typename S>
  double milliseconds() const {
    const auto it = m_milliseconds.find(key<S>());
    if (it != m_milliseconds.end()) {
      return it->second;
    }
    return -1;
  }

protected:
  /**
   * @brief Reset to initial step.
   */
  void reset() {
    m_milliseconds.clear();
  }

private:
  /**
   * @brief Call `get()` on each element of a tuple.
   */
  template <typename STuple, std::size_t... Is>
  void getMultiple(std::index_sequence<Is...>) {
    using mockUnpack = int[];
    (void)mockUnpack {0, (get<std::tuple_element_t<Is, STuple>>(), void(), 0)...};
    // TODO could be done in threads!
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
   * @brief Run step `S` if not done and return its output.
   */
  template <typename S>
  typename S::Return evaluateGet() {
    if (not evaluated<S>()) { // FIXME not thread-safe
      const auto start = std::chrono::high_resolution_clock::now();
      Accessor<S>::evaluate(derived());
      const auto stop = std::chrono::high_resolution_clock::now();
      m_milliseconds[key<S>()] = std::chrono::duration<double, std::milli>(stop - start).count();
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
  template <typename S>
  std::type_index key() const {
    return std::type_index(typeid(S));
  }
  /**
   * @brief The set of performed steps and durations.
   */
  std::map<std::type_index, double> m_milliseconds;
};

} // namespace Framework
} // namespace Phi

#endif // _PHIBOX_STEPPERPIPELINE_H
