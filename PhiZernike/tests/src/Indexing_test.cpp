/**
 * @copyright (C) 2012-2020 Euclid Science Ground Segment
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 3.0 of the License, or (at your option)
 * any later version.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 *
 */

#include "PhiZernike/Indexing.h"

#include <boost/test/unit_test.hpp>

using namespace Phi::Zernike;

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Indexing_test)

//-----------------------------------------------------------------------------

bool operator==(const ZernikeDegrees& lhs, const ZernikeDegrees& rhs) {
  return lhs.radial == rhs.radial && lhs.azimuthal == rhs.azimuthal;
}

BOOST_AUTO_TEST_CASE(conversion_test) {

  struct Indices {
    long ansi;
    long noll;
    long wyant;
    long fringe;
    long radial, azimuthal;
  };

  // From: https://en.wikipedia.org/wiki/Zernike_polynomials
  std::vector<Indices> mapping = {
      {0, 1, 0, 1, 0, 0},
      {1, 3, 2, 3, 1, -1},
      {2, 2, 1, 2, 1, 1},
      {3, 5, 5, 6, 2, -2},
      {4, 4, 3, 4, 2, 0},
      {5, 6, 4, 5, 2, 2},
      {6, 9, 10, 11, 3, -3},
      {7, 7, 7, 8, 3, -1},
      {8, 8, 6, 7, 3, 1},
      {9, 10, 9, 10, 3, 3},
      {10, 15, 17, 18, 4, -4},
      {11, 13, 12, 13, 4, -2},
      {12, 11, 8, 9, 4, 0},
      {13, 12, 11, 12, 4, 2},
      {14, 14, 16, 17, 4, 4}};

  for (const auto& i : mapping) {
    const ZernikeDegrees degrees {i.radial, i.azimuthal};
    const AnsiIndex ansi {i.ansi};
    const NollIndex noll {i.noll};
    const WyantIndex wyant {i.wyant};
    const FringeIndex fringe {i.fringe};

    BOOST_TEST((as<ZernikeDegrees>(degrees) == degrees));

    BOOST_TEST(as<AnsiIndex>(degrees).index == ansi.index);
    BOOST_TEST(as<NollIndex>(degrees).index == noll.index);
    BOOST_TEST(as<WyantIndex>(degrees).index == wyant.index);
    BOOST_TEST(as<FringeIndex>(degrees).index == fringe.index);

    BOOST_TEST((as<ZernikeDegrees>(ansi) == degrees));
    BOOST_TEST((as<ZernikeDegrees>(noll) == degrees));
    BOOST_TEST((as<ZernikeDegrees>(wyant) == degrees));
    BOOST_TEST((as<ZernikeDegrees>(fringe) == degrees));

    BOOST_TEST((as<ZernikeDegrees>(ansi) == degrees));
    BOOST_TEST(as<NollIndex>(ansi).index == noll.index);
    BOOST_TEST(as<WyantIndex>(ansi).index == wyant.index);
    BOOST_TEST(as<FringeIndex>(ansi).index == fringe.index);

    BOOST_TEST(as<AnsiIndex>(ansi).index == ansi.index);
    BOOST_TEST(as<AnsiIndex>(noll).index == ansi.index);
    BOOST_TEST(as<AnsiIndex>(wyant).index == ansi.index);
    BOOST_TEST(as<AnsiIndex>(fringe).index == ansi.index);
  }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
