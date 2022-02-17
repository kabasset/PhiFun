// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef _PHIFOURIER_DFTMEMORY_H
#define _PHIFOURIER_DFTMEMORY_H

#include <complex>
#include <fftw3.h>

namespace Phi {
namespace Fourier {

/**
 * @brief Thread-safe singleton class to ensure proper FFTW memory management.
 * @details
 * This is a Meyer's singleton.
 * The destructor, which is executed once (at the end of the program), calls `fftw_cleanup()`.
 */
class FftwAllocator {
private:
  /**
   * @brief Private constructor.
   */
  FftwAllocator() {}

  /**
   * @brief Destructor which frees.
   */
  ~FftwAllocator() {
    fftw_cleanup();
  }

  /**
   * @brief Get the instance.
   * @details
   * Get the singleton if it exists already or instantiate it otherwise,
   * which triggers cleanup at destruction, i.e. when program ends.
   */
  static FftwAllocator& instantiate() {
    static FftwAllocator allocator;
    return allocator;
  }

public:
  /**
   * @brief Deleted copy constructor.
   */
  FftwAllocator(const FftwAllocator&) = delete;

  /**
   * @brief Deleted copy assignment operator.
   */
  FftwAllocator& operator=(const FftwAllocator&) = delete;

  /**
   * @brief Allocate a buffer's data.
   * @details
   * If `data` is not null, then no allocation is performed and the returned buffer points to `data`.
   */
  template <typename T>
  static T* allocateBuffer(const std::size_t& size, T* data = nullptr) {
    instantiate(); // FIXME is it necessary for buffers?
    return data ? data : (T*)fftw_malloc(sizeof(T) * size);
  }

  /**
   * @brief Free a buffer's data.
   */
  template <typename T>
  static void freeBuffer(T* data) {
    fftw_free(data);
  }

  /**
   * @brief Create a plan.
   * @warning
   * `in` and `out` are filled with garbage.
   */
  template <typename TType, typename TIn, typename TOut>
  static fftw_plan createPlan(TIn& in, TOut& out) {
    instantiate();
    return allocateFftwPlan(in, out);
  }

  /**
   * @brief Destroy a plan.
   */
  static void destroyPlan(fftw_plan plan) {
    fftw_destroy_plan(plan);
  }
};

/// @cond
namespace Internal {

/**
 * @brief Allocate a plan.
 * @warning
 * This does not instantiate the singleton, as opposed to `FftwAllocator::createPlan()`.
 */
template <typename TType, typename TIn, typename TOut>
fftw_plan allocateFftwPlan(TIn& in, TOut& out);

} // namespace Internal
/// @endcond

} // namespace Fourier
} // namespace Phi

#endif // _PHIFOURIER_DFTMEMORY_H
