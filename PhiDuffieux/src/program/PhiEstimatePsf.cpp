// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "EleFits/SifFile.h"
#include "EleFitsData/TestRaster.h"
#include "EleFitsUtils/ProgramOptions.h"
#include "EleFitsValidation/Chronometer.h"
#include "ElementsKernel/ProgramHeaders.h"
#include "PhiDuffieux/MonochromaticSystem.h"
#include "PhiFourier/Dft.h"
#include "PhiZernike/Zernike.h"

#include <map>
#include <string>

using namespace Phi;

using boost::program_options::value;

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

/**
 * @brief Save as a SIF file.
 * @details
 * Does nothing if filename is empty.
 */
template <typename TRaster>
void saveSif(const TRaster& raster, const std::string& filename) {
  if (filename == "") {
    return;
  }
  Euclid::Fits::SifFile f(filename, Euclid::Fits::FileMode::Overwrite);
  f.writeRaster(raster);
  logger.info() << "  See: " << filename;
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

    Euclid::Fits::ProgramOptions options(
        "Compute a monochromatic PSF from a pupil mask and randmom Zernike coefficients.");

    options.named("side", value<long>()->default_value(1024), "Pupil mask side");
    options.named("radius", value<long>()->default_value(256), "Pupil radius");
    options.named("alphas", value<long>()->default_value(40), "Number of Zernike indices");

    options.named("mask", value<std::string>()->default_value("/tmp/mask.fits"), "Pupil mask file");
    options.named("zernike", value<std::string>()->default_value("/tmp/zernike.fits"), "Zernike polynomials file");
    options.named("optics", value<std::string>()->default_value("/tmp/optics.fits"), "Optical PSF file");
    options.named("system", value<std::string>()->default_value("/tmp/system.fits"), "System PSF file");

    return options.asPair();
  }

  /**
   * @brief Run!
   */
  ExitCode mainMethod(std::map<std::string, VariableValue>& args) override {

    const auto maskSide = args["side"].as<long>();
    const auto pupilRadius = args["radius"].as<long>();
    const auto alphaCount = args["alphas"].as<long>();
    const auto maskFilename = args["mask"].as<std::string>();
    const auto zernikeFilename = args["zernike"].as<std::string>();
    const auto opticsFilename = args["optics"].as<std::string>();
    const auto systemFilename = args["system"].as<std::string>();

    using Chrono = Euclid::Fits::Validation::Chronometer<std::chrono::milliseconds>;
    Chrono chrono;

    logger.info("Generating pupil mask...");
    // FIXME load if exists, generate and save otherwise
    chrono.start();
    auto pupil = generatePupil(maskSide, pupilRadius);
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";
    saveSif(pupil, maskFilename);

    logger.info("Generating Zernike polynomials...");
    chrono.start();
    auto zernike = Zernike::basis(maskSide, alphaCount);
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";
    Euclid::Fits::VecRaster<double, 3> zernikeDisp({maskSide, maskSide, alphaCount});
    for (const auto& p : zernikeDisp.domain()) {
      zernikeDisp[p] = zernike[{p[2], p[0], p[1]}];
    }
    saveSif(zernikeDisp, zernikeFilename);

    logger.info("Generating Zernike coefficients...");
    chrono.start();
    std::vector<double> alphas;
    Euclid::Fits::Test::RandomRaster<double, 1>({alphaCount}, -1, 1).moveTo(alphas);
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";
    logger.debug() << "  Coefficients:";
    for (const auto& a : alphas) {
      logger.debug() << "    " << a;
    }

    logger.info("Planning optical DFT and allocating memory...");
    chrono.start();
    Duffieux::MonochromaticOptics optics(.500, maskSide, alphas);
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";

    logger.info("Planning system DFT and allocating memory...");
    chrono.start();
    Duffieux::MonochromaticSystem system(optics);
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";

    chrono.start();
    logger.info("Computing pupil amplitude (complex exp)...");
    optics.evalPupilAmplitude(pupil, zernike);
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";

    logger.info("Computing PSF amplitude (complex DFT)...");
    chrono.start();
    optics.evalPsfAmplitude();
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";

    logger.info("Computing PSF intensity (norm)...");
    chrono.start();
    auto& intensity = optics.evalPsfIntensity();
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";
    saveSif(fftShift(intensity), opticsFilename);

    logger.info("Generating non-optical PSF...");
    chrono.start();
    const Fourier::Position halfShape {(intensity.shape()[0] + 1) / 2, intensity.shape()[1]};
    Fourier::ComplexDftBuffer nonOpticalTf(halfShape);
    // FIXME fill
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";

    logger.info("Computing optical transfer function (real DFT)...");
    chrono.start();
    system.evalOpticalTf();
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";

    logger.info("Computing system transfer function (pixel-wise multiplication)...");
    chrono.start();
    system.evalSystemTf(nonOpticalTf);
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";

    logger.info("Computing system PSF (inverse real DFT)...");
    chrono.start();
    auto& psf = system.evalSystemPsf();
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";
    saveSif(psf, systemFilename); // FIXME fftShift()

    logger.info("Done.");
    return ExitCode::OK;
  }
};

MAIN_FOR(PhiEstimatePsf)
