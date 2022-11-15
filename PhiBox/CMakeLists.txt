CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

elements_subdir(PhiBox)

elements_depends_on_subdirs(LitlCore)
elements_depends_on_subdirs(ElementsKernel)

elements_add_library(PhiBox src/lib/*.cpp
                     INCLUDE_DIRS LitlCore ElementsKernel
                     LINK_LIBRARIES LitlCore ElementsKernel
                     PUBLIC_HEADERS PhiBox)

elements_add_unit_test(ImageProcessing tests/src/ImageProcessing_test.cpp 
                     EXECUTABLE PhiBox_ImageProcessing_test
                     LINK_LIBRARIES PhiBox
                     TYPE Boost)
elements_add_unit_test(SplineIntegrator tests/src/SplineIntegrator_test.cpp 
                     EXECUTABLE PhiBox_SplineIntegrator_test
                     LINK_LIBRARIES PhiBox
                     TYPE Boost)
elements_add_unit_test(StepperPipeline tests/src/StepperPipeline_test.cpp 
                     EXECUTABLE PhiBox_StepperPipeline_test
                     LINK_LIBRARIES PhiBox
                     TYPE Boost)

elements_add_executable(PhiTimeConvolution src/program/PhiTimeConvolution.cpp
                     INCLUDE_DIRS PhiBox
                     LINK_LIBRARIES PhiBox)
elements_add_executable(PhiTimeIntegration src/program/PhiTimeIntegration.cpp
                     INCLUDE_DIRS PhiBox
                     LINK_LIBRARIES PhiBox)

elements_install_conf_files()
