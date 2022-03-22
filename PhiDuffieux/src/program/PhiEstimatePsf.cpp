// Copyright (C) 2022, CNES
// This file is part of PhiFun <github.com/kabasset/PhiFun>
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "EleFits/MefFile.h"
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

    options.named<long>("alphas", "Number of Zernike indices", 40);
    options.named<long>("mask", "Pupil mask diameter", 1024);
    options.named<long>("pupil", "Pupil diameter", 512);
    options.named<long>("psf", "Output PSF diameter", 300);

    options.named<std::string>("output", "Output file", "/tmp/psf.fits");

    return options.asPair();
  }

  /**
   * @brief Run!
   */
  ExitCode mainMethod(std::map<std::string, VariableValue>& args) override {

    const auto alphaCount = args["alphas"].as<long>();
    const auto maskSide = args["mask"].as<long>();
    const auto pupilDiameter = args["pupil"].as<long>();
    const auto psfSide = args["psf"].as<long>();

    Euclid::Fits::MefFile f(args["output"].as<std::string>(), Euclid::Fits::FileMode::Overwrite);

    using Chrono = Euclid::Fits::Validation::Chronometer<std::chrono::milliseconds>;
    Chrono chrono;

    logger.info("Generating pupil mask...");
    // FIXME load if exists, generate and save otherwise
    chrono.start();
    auto pupil = generatePupil(maskSide, pupilDiameter / 2);
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";
    f.appendImage("Pupil mask", {}, pupil);

    logger.info("Generating Zernike polynomials...");
    chrono.start();
    const auto& zernike = Zernike::basis(pupilDiameter / 2, maskSide, alphaCount);
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";
    // Euclid::Fits::VecRaster<double, 3> zernikeDisp({maskSide, maskSide, alphaCount});
    // for (const auto& p : zernikeDisp.domain()) {
    //   zernikeDisp[p] = zernike[{p[2], p[0], p[1]}];
    // }
    // f.appendImage("Zernike basis", {}, zernikeDisp);

    logger.info("Generating Zernike coefficients...");
    chrono.start();
    std::vector<double> alphas;
    Euclid::Fits::Test::RandomRaster<double, 1>({alphaCount}, -1, 1).moveTo(alphas);
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";
    for (std::size_t i = 0; i < alphas.size(); ++i) {
      f.primary().header().write(
          "ALPHA_" + std::to_string(i),
          alphas[i],
          "um",
          "Zernike coefficient " + std::to_string(i));
    }

    logger.info("Planning optical DFT and allocating memory...");
    chrono.start();
    Duffieux::MonochromaticOptics optics(.500, pupil, zernike, alphas);
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";

    logger.info("Planning system DFT and allocating memory...");
    chrono.start();
    Duffieux::MonochromaticSystem system(optics, psfSide);
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";

    logger.info("Computing pupil amplitude (complex exp)...");
    chrono.start();
    const auto& pupilAmplitude = optics.get<Duffieux::PupilAmplitude>();
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";
    f.appendImage("Pupil intensity", {}, norm2(pupilAmplitude));

    logger.info("Computing PSF amplitude (complex DFT)...");
    chrono.start();
    const auto& psfAmplitude = optics.get<Duffieux::PsfAmplitude>();
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";

    logger.info("Computing PSF intensity (norm)...");
    chrono.start();
    const auto& intensity = optics.get<Duffieux::PsfIntensity>();
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";
    f.appendImage("Optical PSF intensity", {}, intensity);

    logger.info("Computing optical transfer function (real DFT)...");
    chrono.start();
    const auto& opticalTf = system.evalOpticalTf();
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";
    f.appendImage("Optical TF intensity", {}, norm2(opticalTf));

    logger.info("Generating non-optical transfer function...");
    chrono.start();
    const Fourier::Position halfShape {(intensity.shape()[0] + 1) / 2, intensity.shape()[1]};
    Fourier::ComplexDftBuffer nonOpticalTf(halfShape);
    std::fill(nonOpticalTf.begin(), nonOpticalTf.end(), 1.); // FIXME
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";
    f.appendImage("Non-optical TF intensity", {}, norm2(nonOpticalTf));

    logger.info("Computing system transfer function (pixel-wise multiplication)...");
    chrono.start();
    const auto& stf = system.evalSystemTf(nonOpticalTf);
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";
    f.appendImage("System TF intensity", {}, norm2(stf));

    logger.info("Wrapping system transfer function (warping)...");
    chrono.start();
    const auto& warpedTf = system.warpSystemTf(1, 1, 1, 1);
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";
    f.appendImage("Wrapped system TF intensity", {}, norm2(warpedTf));

    logger.info("Computing system PSF (inverse real DFT)...");
    chrono.start();
    const auto& psf = system.evalSystemPsf();
    chrono.stop();
    logger.info() << "  " << chrono.last().count() << "ms";
    f.appendImage("System PSF", {}, psf);

    logger.info("Done.");
    return ExitCode::OK;
  }
};

MAIN_FOR(PhiEstimatePsf)
