// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "EleFitsData/Raster.h"
#include "EleFitsData/TestRaster.h"
#include "EleFitsUtils/ProgramOptions.h"
#include "EleFitsValidation/Chronometer.h"
#include "ElementsKernel/ProgramHeaders.h"
#include "PhiBox/ImageProcessing.h"
#include "PhiBox/SplineIntegrator.h"

#include <map>
#include <string>

using namespace Phi;

using boost::program_options::value; // FIXME rm

using Raster2D = Euclid::Fits::VecRaster<std::complex<double>, 2>;
using Raster3D = Euclid::Fits::VecRaster<std::complex<double>, 3>;
using View2D = Euclid::Fits::PtrRaster<std::complex<double>, 2>;
using Chrono = Euclid::Fits::Validation::Chronometer<std::chrono::milliseconds>;

void warp(const Raster3D& input, Raster3D& output) {
  const auto width = output.shape()[0];
  const auto height = output.shape()[1];
  const auto depth = output.shape()[2];
  const double xFactor = double(input.shape()[0] - 1) / (width - 1);
  const double yFactor = double(input.shape()[1] - 1) / (height - 1);
  auto* it = output.begin();
  for (long l = 0; l < depth; ++l) {
    for (long j = 0; j < height; ++j) {
      const double v = yFactor * j;
      for (long i = 0; i < width; ++i, ++it) {
        const double u = xFactor * i;
        *it = Image2D::bilinear<std::complex<double>>(input, u / l, v / l, l);
      }
    }
  }
}

std::vector<double> linspace(double min, double sup, long count) {
  std::vector<double> x(count);
  for (long i = 0; i < count; ++i) {
    x[i] = min + i * (sup - min) / count;
  }
  return x;
}

Raster2D
integrate(const Raster3D& input, const Spline::SplineIntegrator& integrator, const std::vector<double>& weights) {
  const auto width = input.shape()[0];
  const auto height = input.shape()[1];
  const auto size = width * height;
  Raster2D output({width, height});
  const auto depth = input.shape()[2];
  std::vector<std::complex<double>> y(depth);
  std::vector<std::complex<double>> z(depth);
  auto oIt = output.begin();
  auto iIt = input.begin();
  const auto d = iIt - oIt;
  for (long ij = 0; ij < size; ++ij, ++oIt, iIt = oIt + d) {
    auto yIt = y.data();
    for (long l = 0; l < depth; ++l, ++yIt, iIt += size) {
      *yIt = *iIt;
    }
    z = integrator.knotZ(y.data());
    *oIt = integrator.integrate(y.data(), z.data(), weights.data());
  }
  return output;
}

void warpIntegrate(
    const Raster3D& input,
    Raster2D& output,
    const Spline::SplineIntegrator& integrator,
    const std::vector<double>& weights) {
  const auto width = output.shape()[0];
  const auto height = output.shape()[1];
  const double xFactor = double(input.shape()[0] - 1) / (width - 1);
  const double yFactor = double(input.shape()[1] - 1) / (height - 1);
  const auto depth = input.shape()[2];
  std::vector<std::complex<double>> y(depth);
  std::vector<std::complex<double>> z(depth);
  auto* it = output.begin();
  for (long j = 0; j < height; ++j) {
    const double v = yFactor * j;
    for (long i = 0; i < width; ++i, ++it) {
      const double u = xFactor * i;
      for (long l = 0; l < depth; ++l) {
        y[l] = Image2D::bilinear<std::complex<double>>(input, u / l, v / l, l);
      }
      z = integrator.knotZ(y.data());
      *it = integrator.integrate(y.data(), z.data(), weights.data());
    }
  }
}

class PhiTimeIntegration : public Elements::Program {

public:
  std::pair<OptionsDescription, PositionalOptionsDescription> defineProgramArguments() override {
    Euclid::Fits::ProgramOptions options("Wrap and integrate a set of random monochromatic TFs into a broadband TF.");
    options.named("input", "Monochromatic TF diameter", 1024);
    options.named("output", "Broadband TF diameter", 512);
    options.named("lambdas", "Number of wavelengths", 40);
    options.named("steps", "Number of integration steps", 200);
    options.flag("combine", "Combine the interpolation and integration steps");
    return options.asPair();
  }

  ExitCode mainMethod(std::map<std::string, VariableValue>& args) override {
    Logging logger = Logging::getLogger("PhiTimeIntegration");

    const auto inputSide = args["input"].as<long>();
    const auto outputSide = args["output"].as<long>();
    const auto lambdaCount = args["lambdas"].as<long>();
    const auto stepCount = args["steps"].as<long>();

    Chrono chrono;

    logger.info() << "Generating values...";
    chrono.start();
    Euclid::Fits::Test::RandomRaster<std::complex<double>, 3> input({inputSide / 2 + 1, inputSide, lambdaCount});
    chrono.stop();
    logger.info() << "  " << inputSide / 2 + 1 << "x" << inputSide << "x" << lambdaCount << "px";
    logger.info() << "  " << chrono.last().count() << "ms";

    logger.info() << "Pre-computing integration coefficients...";
    chrono.start();
    Spline::SplineIntegrator integrator(linspace(500., 900., lambdaCount), linspace(500., 900., stepCount));
    const std::vector<double> weights(integrator.interpolationX().size(), 1.);
    chrono.stop();
    logger.info()
        << "  " << integrator.knotX().size() << " to " << integrator.interpolationX().size() << " values per pixel";
    logger.info() << "  " << chrono.last().count() << "ms";

    if (not args["combine"].as<bool>()) {
      logger.info() << "Wrapping...";
      chrono.start();
      Raster3D output({outputSide / 2 + 1, outputSide, lambdaCount});
      warp(input, output);
      chrono.stop();
      logger.info() << "  " << outputSide / 2 + 1 << "x" << outputSide << "x" << lambdaCount << "px";
      logger.info() << "  " << chrono.last().count() << "ms";

      logger.info() << "Integrating...";
      chrono.start();
      auto broadband = integrate(output, integrator, weights);
      chrono.stop();
      logger.info() << "  " << chrono.last().count() << "ms";
    } else {
      logger.info() << "Wrapping and integrating...";
      chrono.start();
      Raster2D broadband({outputSide / 2 + 1, outputSide});
      warpIntegrate(input, broadband, integrator, weights);
      chrono.stop();
      logger.info() << "  " << chrono.last().count() << "ms";
    }

    return ExitCode::OK;
  }
};

MAIN_FOR(PhiTimeIntegration)
