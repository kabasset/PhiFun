// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef _PHIFOURIER_DFTBUFFER_H
#define _PHIFOURIER_DFTBUFFER_H

#include "LitlCore/Raster.h"
#include "PhiFourier/DftMemory.h"

namespace Phi {
namespace Fourier {

/**
 * @brief Data holder for FFTW's allocated buffers.
 * @details
 * Data can be either owned, or shared and owned by another object.
 */
template <typename T>
struct DftBufferData {

public:
  /**
   * @brief The concrete container type.
   */
  using Container = T*;

  /**
   * @brief Constructor.
   * @details
   * Allocate some memory if `data = nullptr`.
   */
  DftBufferData(std::size_t size, T* data = nullptr) :
      m_shared(data), m_size(size), m_container(data ? data : FftwAllocator::allocateBuffer<T>(m_size)) {}

  /**
   * @brief Non-copyable.
   */
  DftBufferData(const DftBufferData&) = delete;

  /**
   * @brief Movable.
   */
  DftBufferData(DftBufferData&& other) :
      m_shared(other.m_shared), m_size(other.m_size), m_container(other.m_container) {
    other.m_shared = false;
    other.m_size = 0;
    other.m_container = nullptr;
  }

  /**
   * @brief Non-copyable.
   */
  DftBufferData& operator=(const DftBufferData&) = delete;

  /**
   * @brief Movable.
   */
  DftBufferData& operator=(DftBufferData&& other) {
    if (this != &other) {
      m_shared = other.m_shared;
      m_size = other.m_size;
      m_container = other.m_container;
      other.m_shared = false;
      other.m_size = 0;
      other.m_container = nullptr;
    }
    return *this;
  }

  /**
   * @brief Destructor.
   * @details
   * Free memory if needed.
   */
  ~DftBufferData() {
    if (owns()) {
      FftwAllocator::freeBuffer(m_container);
      m_container = nullptr;
    }
  }

  /**
   * @brief Get the data size.
   */
  std::size_t size() const {
    return m_size;
  }

  /**
   * @brief Get the data pointer.
   */
  const T* data() const {
    return m_container;
  }

  /**
   * @brief Check whether the data is owned by this object.
   */
  bool owns() const {
    return m_container && not m_shared;
  }

protected:
  /**
   * @brief Is the data shared?
   */
  bool m_shared;

  /**
   * @brief The data size.
   */
  std::size_t m_size;

  /**
   * @brief The data pointer.
   */
  T* m_container;
};

/**
 * @brief A 2D position.
 */
using Position = Litl::Position<2>;

/**
 * @brief Input or output buffer of a `DftPlan`.
 */
template <typename T>
using DftBuffer = Litl::Raster<T, 2, DftBufferData<T>>;

/**
 * @brief Specialization for `double`.
 */
using RealDftBuffer = DftBuffer<double>;

/**
 * @brief Specialization for `std::complex<double>`.
 */
using ComplexDftBuffer = DftBuffer<std::complex<double>>;

/**
 * @brief Multiply by a buffer.
 */
template <typename T, typename TRaster>
DftBuffer<T>& operator*=(DftBuffer<T>& lhs, const TRaster& rhs) {
  std::transform(lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(), [](const T& l, const T& r) {
    return l * r;
  });
  return lhs;
}

} // namespace Fourier
} // namespace Phi

#endif // _PHIFOURIER_DFTBUFFER_H
