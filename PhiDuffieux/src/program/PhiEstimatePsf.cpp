// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "EleFits/MefFile.h"
#include "EleFitsData/TestRaster.h"
#include "EleFitsUtils/ProgramOptions.h"
#include "EleFitsValidation/Chronometer.h"
#include "ElementsKernel/ProgramHeaders.h"
#include "PhiBox/SplineIntegrator.h"
#include "PhiDuffieux/MonochromaticSystem.h"
#include "PhiFourier/Dft.h"
#include "PhiZernike/Zernike.h"

#include <map>
#include <string>

using namespace Phi;

using boost::program_options::value; // FIXME rm

static Elements::Logging logger = Elements::Logging::getLogger("PhiEstimatePsf");

/**
 * @brief Generate a circular pupil mask.
 */
Euclid::Fits::VecRaster<double> generatePupil(long maskSide, long pupilRadius) {
  Euclid::Fits::VecRaster<double> pupil({maskSide, maskSide});
  const auto maskRadius = maskSide / 2;
  const auto pupilRadiusSquared = pupilRadius * pupilRadius;
  for (const auto& p : pupil.domain()) {
    const auto u = double(p[0] - maskRadius);
    const auto v = double(p[1] - maskRadius);
    if (std::abs(u * u + v * v) < pupilRadiusSquared) {
      pupil[p] = 1.;
    }
  }
  return pupil;
}

std::string wavelengthString(double lambda) {
  return std::to_string(int(lambda + .5)) + " nm";
}

/**
 * @brief The program.
 */
class PhiEstimatePsf : public Elements::Program {

public:
  /**
   * @brief The options.
   */
  std::pair<OptionsDescription, PositionalOptionsDescription> defineProgramArguments() override {

    Euclid::Fits::ProgramOptions options("Compute a broadband PSF from a pupil mask and randmom Zernike coefficients.");

    options.named("alphas", "Number of Zernike indices", 45L);
    options.named("lambdas", "Number of wavelengths", 16L);
    options.named("mask", "Pupil mask diameter", 1024L);
    options.named("pupil", "Pupil diameter", 512L);
    options.named("psf", "Output PSF diameter", 300L);

    options.named("output", "Output file", std::string("/tmp/broadband.fits"));

    return options.asPair();
  }

  /**
   * @brief Run!
   */
  ExitCode mainMethod(std::map<std::string, VariableValue>& args) override {

    const auto alphaCount = args["alphas"].as<long>();
    const auto lambdaCount = args["lambdas"].as<long>();
    const auto maskSide = args["mask"].as<long>();
    const auto pupilDiameter = args["pupil"].as<long>();
    const auto psfSide = args["psf"].as<long>();

    Euclid::Fits::MefFile f(args["output"].as<std::string>(), Euclid::Fits::FileMode::Overwrite);

    using Chrono = Euclid::Fits::Validation::Chronometer<std::chrono::milliseconds>;
    Chrono chrono;

    logger.info("Initialization...");
    const auto pupil = generatePupil(maskSide, pupilDiameter / 2);
    logger.info("  Pupil mask done.");
    const auto zernike = Zernike::basis(pupilDiameter / 2, maskSide, alphaCount);
    logger.info("  Zernike basis done.");
    std::vector<double> alphas;
    Euclid::Fits::Test::RandomRaster<double, 1>({alphaCount}, -1000, 1000).moveTo(alphas);
    logger.info("  Zernike coefficients done.");
    Fourier::ComplexDftBuffer nonOpticalTf({maskSide / 2 + 1, maskSide});
    std::fill(nonOpticalTf.begin(), nonOpticalTf.end(), 1.); // FIXME
    logger.info("  Non-optical transfer function done.");

    logger.info("Planning DFTs and allocating buffers...");
    std::vector<Duffieux::MonochromaticSystem> systems;
    systems.reserve(lambdaCount);
    for (double lambda : Spline::linspace(500., 900., lambdaCount)) {
      chrono.start();
      systems.emplace_back(
          Duffieux::MonochromaticOptics::Params {lambda, pupil, zernike, alphas},
          Duffieux::MonochromaticSystem::Params {{psfSide, psfSide}, nonOpticalTf, {.0001, 0, 0, .0001}});
      chrono.stop();
      logger.info() << "  " << wavelengthString(lambda) << " system done in: " << chrono.last().count() << " ms";
    }

    logger.info("Monochromatic PSF computation...");
    chrono.start();
#pragma omp parallel for
    for (auto& system : systems) {
      system.get<Duffieux::WarpedSystemPsf>();
      logger.info()
          << "  " << wavelengthString(system.wavelength()) << " nm PSF computation: " << system.milliseconds() << " ms";
    }
    chrono.stop();
    logger.info() << "Overall computation: " << chrono.last().count() << " ms";

    return ExitCode::OK;
  }
};

MAIN_FOR(PhiEstimatePsf)
