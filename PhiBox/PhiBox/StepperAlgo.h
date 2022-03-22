// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef _PHIBOX_STEPPERALGO_H
#define _PHIBOX_STEPPERALGO_H

#include <set>
#include <typeindex>

namespace Phi {
namespace Framework {

/**
 * @brief Boolean wrapper class which test whether a given step has a prerequisite.
 */
template <typename S>
struct HasPrerequisite {
  static constexpr bool value = not std::is_same<typename S::Prerequisite, void>::value;
};

/**
 * @brief An algorithm which can be run step-by-step using lazy evaluation.
 * @details
 * Only one method is provided, `get<S>()`, which returns the value of step `S`.
 * 
 * This class relies on the CRTP.
 * Child classes should provide a specialization of the following methods for each step `S`:
 * - `void doEvaluate<S>()`, which evaluates `S` assuming prerequisites were already computed;
 * - `S::Return doGet<S>()`, which returns the computed value of `S`.
 * 
 * A step `S` is a class which contains the following types:
 * - `Return` is the value type of the output of step `S`;
 * - `Prerequisite` is the step which must be run prior to `S`, or `void` if there is no prerequisite.
 */
template <typename Algo>
class StepperAlgo {

public:
  /**
   * @brief Lazy evaluation of step `S` which triggers prerequisite evaluation if needed.
   */
  template <typename S>
  std::enable_if_t<HasPrerequisite<S>::value, typename S::Return> get() {
    get<typename S::Prerequisite>();
    return evaluateGet<S>();
  }

  /**
   * @brief Evaluation of step `S` which has no prerequisite.
   */
  template <typename S>
  std::enable_if_t<not HasPrerequisite<S>::value, typename S::Return> get() {
    return evaluateGet<S>();
  }

  /**
   * @brief Reset to initial step.
   */
  void reset() {
    m_done.clear();
  }

private:
  /**
   * @brief Access to protected methods of `Algo`.
   */
  template <typename S>
  struct Accessor : Algo {

    /**
     * @brief Call `algo.doEvaluate<S>()`.
     */
    static void evaluate(Algo& algo) {
      void (Algo::*fn)() = &Accessor::template doEvaluate<S>;
      (algo.*fn)();
    }

    /**
     * @brief Call `algo.doGet<S>()`.
     */
    static typename S::Return get(Algo& algo) {
      typename S::Return (Algo::*fn)() const = &Accessor::template doGet<S>;
      return (algo.*fn)();
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
   * @brief Cast as `Algo`.
   */
  inline Algo& derived() {
    return static_cast<Algo&>(*this);
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
