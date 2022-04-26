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

class PhiTimeConvolution : public Elements::Program {

public:
  std::pair<OptionsDescription, PositionalOptionsDescription> defineProgramArguments() override {
    Euclid::Fits::ProgramOptions options("Convolve by a separable filter with early or lazy decimation.");
    options.named("input", "Input map diameter", 1024);
    options.named("kernel", "Kernel size", 24);
    return options.asPair();
  }

  ExitCode mainMethod(std::map<std::string, VariableValue>& args) override {

    Logging logger = Logging::getLogger("PhiTimeConvolution");

    return ExitCode::OK;
  }
};

MAIN_FOR(PhiTimeConvolution)
