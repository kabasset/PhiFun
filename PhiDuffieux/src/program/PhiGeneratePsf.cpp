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

static Elements::Logging logger = Elements::Logging::getLogger("PhiGeneratePsf");

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

/**
 * @brief The program.
 */
class PhiGeneratePsf : public Elements::Program {

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

    logger.info("Generating pupil mask...");
    const auto pupil = generatePupil(maskSide, pupilDiameter / 2);
    const auto zernike = Zernike::basis(pupilDiameter / 2, maskSide, alphaCount);

    logger.info("Generating Zernike coefficients...");
    chrono.start();
    std::vector<double> alphas;
    Euclid::Fits::Test::RandomRaster<double, 1>({alphaCount}, -1000, 1000).moveTo(alphas);
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";
    for (std::size_t i = 0; i < alphas.size(); ++i) {
      f.primary().header().write(
          "ALPHA_" + std::to_string(i),
          alphas[i],
          "nm",
          "Zernike coefficient " + std::to_string(i));
    }

    logger.info("Generating non-optical transfer function...");
    Fourier::ComplexDftBuffer nonOpticalTf({maskSide / 2 + 1, maskSide});
    std::fill(nonOpticalTf.begin(), nonOpticalTf.end(), 1.); // FIXME

    logger.info("Planning monochromatic DFTs and allocating memory...");
    chrono.start();
    Duffieux::MonochromaticOptics optics(.500, pupil, zernike, alphas);
    Duffieux::MonochromaticSystem system(optics, psfSide);
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";

    // Loop over lambdas to build the monochromatic system TFs
    for (double lambda : Spline::linspace(500., 900., lambdaCount)) {
      optics.updateLambda(lambda);
      const auto lambdaStr = std::to_string(int(lambda + .5)) + "nm";
      logger.info() << "Lambda = " << lambdaStr;

      logger.info("  Computing PSF intensity (complex exp, complex DFT, norm)...");
      chrono.start();
      const auto& intensity = optics.psfIntensity();
      chrono.stop();
      logger.info() << "    " << chrono.last().count() << "ms";
      f.appendImage(lambdaStr + " optical PSF", {}, intensity);

      logger.info("  Computing optical transfer function (real DFT, multiplication)...");
      chrono.start();
      system.evalOpticalTf();
      const auto& stf = system.evalSystemTf(nonOpticalTf);
      chrono.stop();
      logger.info() << "    " << chrono.last().count() << "ms";
      f.appendImage(lambdaStr + " system TF", {}, norm2(stf));
    }

    logger.info("Wrapping system transfer function (bilinear interpolation)...");
    chrono.start();
    const auto& warpedTf = system.warpSystemTf();
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";
    f.appendImage("Wrapped system TF intensity", {}, norm2(warpedTf));

    logger.info("Computing system PSF (inverse real DFT)...");
    chrono.start();
    const auto& psf = system.evalSystemPsf();
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";
    f.appendImage("System PSF", {}, psf);

    logger.info() << "See " + f.filename();
    return ExitCode::OK;
  }
};

MAIN_FOR(PhiGeneratePsf)
