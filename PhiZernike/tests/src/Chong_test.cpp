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

#include "EleFits/SifFile.h" // FIXME rm
#include "EleFitsValidation/Chronometer.h" // FIXME rm
#include "PhiZernike/Chong.h"

#include <boost/test/unit_test.hpp>

using namespace Phi::Zernike;

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Chong_test)

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(save_zernikes_test) {
  Euclid::Fits::Validation::Chronometer<std::chrono::milliseconds> chrono;
  constexpr long side = 1024;
  constexpr long orders = 10;
  chrono.start();
  const auto zxy = chongBasis(side, orders);
  chrono.stop();
  printf("Elapsed: %lims\n", chrono.last().count());
  Euclid::Fits::VecRaster<double, 3> xyz({side, side, zxy.shape()[0]});
  for (const auto& p : xyz.domain()) {
    xyz[p] = zxy[{p[2], p[1], p[0]}];
  }
  Euclid::Fits::SifFile f("/tmp/zernike.fits", Euclid::Fits::FileMode::Overwrite);
  f.writeRaster(xyz);
}

//-----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
