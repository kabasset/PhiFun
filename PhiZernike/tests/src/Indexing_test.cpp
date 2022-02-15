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

BOOST_AUTO_TEST_CASE(conversion_test) {

  struct Indices {
    long radial, azimuthal;
    long ansi;
    long noll;
  };

  std::vector<Indices> mapping = {
      {0, 0, 0, 1},
      {1, -1, 1, 3},
      {1, 1, 2, 2},
      {2, -2, 3, 5},
      {2, 0, 4, 4},
      {2, 2, 5, 6},
      {3, -3, 6, 9},
      {3, -1, 7, 7},
      {3, 1, 8, 8},
      {3, 3, 9, 10},
      {4, -4, 10, 15},
      {4, -2, 11, 13},
      {4, 0, 12, 11},
      {4, 2, 13, 12},
      {4, 4, 14, 14}};

  for (const auto& i : mapping) {
    BOOST_TEST(as<AnsiIndex>(ZernikeDegrees {i.radial, i.azimuthal}).index == i.ansi);
    BOOST_TEST(as<AnsiIndex>(NollIndex {i.noll}).index == i.ansi);
  }
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
