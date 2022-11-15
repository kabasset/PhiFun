// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "EleFits/MefFile.h"
#include "ElementsKernel/ProgramHeaders.h"
#include "LitlRun/Chronometer.h"
#include "LitlRun/ProgramOptions.h"
#include "PhiBox/SplineIntegrator.h"
#include "PhiDuffieux/MonochromaticSystem.h"
#include "PhiFourier/Dft.h"
#include "PhiZernike/Zernike.h"

#include <map>
#include <string>

using namespace Phi;

using boost::program_options::value; // FIXME rm

static Elements::Logging logger = Elements::Logging::getLogger("PhiGeneratePsf");

// FIXME rm when EleFits moves to Litl
template <typename TRaster>
void appendImage(Euclid::Fits::MefFile& f, const std::string& name, const TRaster& raster) {
  f.appendImage(name, {}, Euclid::Fits::makeRaster(raster.data(), raster.shape()[0], raster.shape()[1]));
}

/**
 * @brief Generate a circular pupil mask.
 */
Litl::Raster<double> generatePupil(long maskSide, long pupilRadius) {
  Litl::Raster<double> pupil({maskSide, maskSide});
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

    Litl::ProgramOptions options("Compute a broadband PSF from a pupil mask and randmom Zernike coefficients.");

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

    using Chrono = Litl::Chronometer<std::chrono::milliseconds>;
    Chrono chrono;

    logger.info("Generating pupil mask...");
    chrono.start();
    const auto pupil = generatePupil(maskSide, pupilDiameter / 2);
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";

    logger.info("Generating Zernike basis...");
    chrono.start();
    const auto zernike = Zernike::basis(pupilDiameter / 2, maskSide, alphaCount);
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";

    logger.info("Generating random Zernike coefficients...");
    chrono.start();
    std::vector<double> alphas;
    Litl::Raster<double, 1>({alphaCount}).generate(Litl::UniformNoise<double>(-1000, 1000)).moveTo(alphas);
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";
    for (std::size_t i = 0; i < alphas.size(); ++i) {
      f.primary().header().write(
          "ALPHA_" + std::to_string(i),
          alphas[i],
          "nm",
          "Zernike coefficient " + std::to_string(i));
    }

    logger.info("Generating identity non-optical transfer function...");
    chrono.start();
    Fourier::ComplexDftBuffer nonOpticalTf({maskSide / 2 + 1, maskSide});
    std::fill(nonOpticalTf.begin(), nonOpticalTf.end(), 1.); // FIXME
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";

    logger.info("Planning monochromatic DFTs and allocating memory...");
    chrono.start();
    Duffieux::MonochromaticSystem system(
        {.500, pupil.shape(), pupil.data(), zernike.data(), alphas},
        {{psfSide, psfSide}, nonOpticalTf.data(), {.0001, 0, 0, .0001}});
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";

    // Loop over lambdas to build the monochromatic system TFs
    for (double lambda : Spline::linspace(500., 900., lambdaCount)) {
      system.update(lambda);
      const auto lambdaStr = std::to_string(int(lambda + .5)) + " nm";
      logger.info() << "Lambda = " << lambdaStr;

      logger.info("  Computing PSF intensity (complex exp, complex DFT, norm)...");
      const auto& intensity = system.get<Duffieux::PsfIntensity>();
      logger.info() << "    Complex exp: " << system.milliseconds<Duffieux::PupilAmplitude>() << " ms";
      logger.info() << "    Complex DFT: " << system.milliseconds<Duffieux::PsfAmplitude>() << " ms";
      logger.info() << "    Norm squared: " << system.milliseconds<Duffieux::PsfIntensity>() << " ms";
      appendImage(f, lambdaStr + " optical PSF", intensity);

      logger.info("  Computing optical transfer function...");
      const auto& stf = system.get<Duffieux::SystemTf>();
      logger.info() << "    Real DFT, complex multiplication: " << system.milliseconds<Duffieux::SystemTf>() << " ms";
      appendImage(f, lambdaStr + " system TF", norm2(stf));

      logger.info("  Wrapping system transfer function...");
      const auto& warpedTf = system.get<Duffieux::WarpedSystemTf>();
      logger.info() << "    Bilinear interpolation: " << system.milliseconds<Duffieux::WarpedSystemTf>() << " ms";
      appendImage(f, lambdaStr + " warped system TF intensity", norm2(warpedTf));

      logger.info("  Computing system PSF...");
      const auto& psf = system.get<Duffieux::WarpedSystemPsf>();
      logger.info() << "    Inverse real DFT: " << system.milliseconds<Duffieux::WarpedSystemPsf>() << " ms";
      appendImage(f, lambdaStr + " system PSF", psf);

      logger.info() << "  Overall: " << system.milliseconds() << " ms";
    }

    logger.info() << "See " + f.filename();
    return ExitCode::OK;
  }
};

MAIN_FOR(PhiGeneratePsf)
