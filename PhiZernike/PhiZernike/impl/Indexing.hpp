// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef _PHIZERNIKE_IMPL_INDEXING_HPP
#define _PHIZERNIKE_IMPL_INDEXING_HPP

namespace Phi {
namespace Zernike {

namespace Internal {

constexpr long constexprBoundedSqrt(long res, long lower, long upper) {
  if (lower == upper) {
    return upper;
  } else {
    const auto mid = (upper + lower) / 2;
    if (mid * mid >= res) {
      return constexprBoundedSqrt(res, lower, mid);
    } else {
      return constexprBoundedSqrt(res, mid + 1, upper);
    }
  }
}

static constexpr long constexprSqrt(long res) {
  return constexprBoundedSqrt(res, 1, res);
}

template <typename TTo, typename TFrom>
struct IndexConverter {
  static constexpr TTo convert(const TFrom& from);
};

template <>
struct IndexConverter<ZernikeDegrees, ZernikeDegrees> {
  static constexpr ZernikeDegrees convert(const ZernikeDegrees& degrees);
};

template <typename TTo>
struct IndexConverter<TTo, ZernikeDegrees> {
  static constexpr TTo convert(const ZernikeDegrees& from);
};

template <typename TFrom>
struct IndexConverter<ZernikeDegrees, TFrom> {
  static constexpr ZernikeDegrees convert(const TFrom& from);
};

// Idendity

constexpr ZernikeDegrees IndexConverter<ZernikeDegrees, ZernikeDegrees>::convert(const ZernikeDegrees& degrees) {
  return degrees;
}

// From Zernike

template <>
constexpr AnsiIndex IndexConverter<AnsiIndex, ZernikeDegrees>::convert(const ZernikeDegrees& degrees) {
  return AnsiIndex {(degrees.radial * (degrees.radial + 2) + degrees.azimuthal) / 2};
}

template <>
constexpr NollIndex IndexConverter<NollIndex, ZernikeDegrees>::convert(const ZernikeDegrees& degrees) {
  const long n = degrees.radial;
  long l = degrees.azimuthal;
  long rhs = 0;
  if (n % 4 <= 1) {
    if (l <= 0) {
      l = -l;
      rhs = 1;
    }
  } else if (l < 0) {
    l = -l;

  } else {
    rhs = 1;
  }
  return NollIndex {n * (n + 1) / 2 + l + rhs};
}

template <>
constexpr FringeIndex IndexConverter<FringeIndex, ZernikeDegrees>::convert(const ZernikeDegrees& degrees) {
  long abs = degrees.azimuthal;
  long sgn = 1;
  if (abs < 0) {
    abs = -abs;
    sgn = -1;
  } else if (abs == 0) {
    sgn = 0;
  }
  const long brace = 1 + (degrees.radial + abs) / 2;
  return FringeIndex {brace * brace - 2 * abs + (1 - sgn) / 2};
}

template <>
constexpr WyantIndex IndexConverter<WyantIndex, ZernikeDegrees>::convert(const ZernikeDegrees& degrees) {
  return WyantIndex {IndexConverter<FringeIndex, ZernikeDegrees>::convert(degrees).index - 1};
}

// To Zernike

template <>
constexpr ZernikeDegrees IndexConverter<ZernikeDegrees, AnsiIndex>::convert(const AnsiIndex& ansi) {
  const auto j = ansi.index;
  long radial = int((constexprSqrt(8 * j + 9) - 2) / 2);
  const long azimuthal = 2 * j - radial * (radial + 2);
  return ZernikeDegrees {radial, azimuthal};
}

template <>
constexpr ZernikeDegrees IndexConverter<ZernikeDegrees, NollIndex>::convert(const NollIndex& noll) {
  long radial = 0;
  long i = noll.index - 1;
  while (i > radial) {
    ++radial;
    i -= radial;
  }
  long azimuthal = (radial % 2) + 2 * long((i + ((radial + 1) % 2)) / 2);
  if (noll.index % 2 == 0) {
    return ZernikeDegrees {radial, azimuthal};
  } else {
    return ZernikeDegrees {radial, -azimuthal};
  }
}

template <>
constexpr ZernikeDegrees IndexConverter<ZernikeDegrees, FringeIndex>::convert(const FringeIndex& fringe) {
  // Brute force!
  for (long n = 0;; ++n) {
    for (long m = -n; m <= n; m += 2) {
      ZernikeDegrees candidate {n, m};
      if (IndexConverter<FringeIndex, ZernikeDegrees>::convert(candidate).index == fringe.index) {
        return candidate;
      }
    }
  }
}

template <>
constexpr ZernikeDegrees IndexConverter<ZernikeDegrees, WyantIndex>::convert(const WyantIndex& wyant) {
  return IndexConverter<ZernikeDegrees, FringeIndex>::convert(FringeIndex {wyant.index + 1});
}

} // namespace Internal

template <typename TTo, typename TFrom>
constexpr TTo as(const TFrom& from) {
  return Internal::IndexConverter<TTo, ZernikeDegrees>::convert(
      Internal::IndexConverter<ZernikeDegrees, TFrom>::convert(from));
}

} // namespace Zernike
} // namespace Phi

#endif // _PHIZERNIKE_IMPL_INDEXING_HPP
