// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef _PHIBOX_SAMPLING_H
#define _PHIBOX_SAMPLING_H

#include "EleFitsData/Region.h"

namespace Phi {

/**
 * @brief Linear sampling.
 */
template <typename T>
class ScalarLinSampling {

public:
  ScalarLinSampling(T front = 0, T back = -1, T step = 1, T stride = 1) :
      m_size(size), m_front(front), m_back(back), m_step(step), m_stride(stride) {}

  parse(const std::string& spec) {
    auto tokens = Euclid::Fits::String::split(spec, ":");
  }

  std::size_t size() const {
    return (m_back - m_front + m_step) / m_step;
  }

  T front() const {
    return m_front;
  }

  ScalarLinSampling& front(T value) {
    m_front = value;
    return *this;
  }

  T back() const {
    return m_back;
  }

  ScalarLinSampling& back(T value) {
    m_back = value;
    return *this;
  }

  T step() const {
    return m_step;
  }

  ScalarLinSampling& step(T value) {
    m_step = value;
    return *this;
  }

  T stride() const {
    return m_stride;
  }

  ScalarLinSampling& stride(T value) {
    m_stride = value;
    return *this;
  }

  std::vector <

      private : T m_size;
  T m_front;
  T m_back;
  T m_step;
  T m_stride;
};

template <typename T>
class ScalarLinSamples {

public:
  using Value = T;

  template <typename TSamples>
  class Iterator : public std::iterator<std::input_iterator_tag, T> {

  public:
    using Value = typename TSamples::Value;

    Iterator(TSamples& sampling, T index) : m_sampling(sampling), m_it(sampling.data() + index * sampling.stride()) {}

    Value& operator*() {
      return *m_it;
    }

    Value* operator->() {
      return m_it;
    }

    Iterator& operator++() {
      return *this += 1;
    }

    Iterator operator++(int) {
      auto res = *this;
      ++res;
    }

    Iterator& operator+=(T n) {
      m_it += m_sampling.step() * m_sampling.stride() * n;
      return *this;
    }

    Iterator& operator-=(T n) {
      *this += -n;
      return *this;
    }

    bool operator==(const Iterator& rhs) const {
      return m_it == rhs.m_it;
    }

    bool operator!=(const Iterator& rhs) const {
      return not(*this == rhs);
    }

    Iterator& operator=(Value* it) {
      m_it = it;
      return *this;
    }

  private:
    TSamples& m_sampling;
    Value* m_it;
  };

  const T* data() const {
    return m_data;
  }

  T* data() {
    return m_data;
  }

  ScalarLinSampling& data(T* d) {
    m_data = d;
    return *this;
  }

  Iterator<const ScalarLinSamples<const T>> begin() const {
    return {*this, m_front};
  }

  Iterator<const ScalarLinSamples<const T>> end() const {
    return {*this, m_back + 1};
  }

  Iterator<ScalarLinSamples<T>> begin() {
    return {*this, m_front};
  }

  Iterator<ScalarLinSamples<T>> end() {
    return {*this, m_back + 1};
  }

private:
  ScalarLinSampling m_sampling;
  T* m_data;
};

} // namespace Phi

#endif // _PHIBOX_SAMPLING_H
