CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

#===============================================================================
# Load elements_subdir macro here
# Examples:
#   For declaring a project module:
#         elements_subdir(ElementsExamples)
#===============================================================================
elements_subdir(PhiBox)

#===============================================================================
# Load elements_depends_on_subdirs macro here
#   For creating a dependency onto an other accessible module
#         elements_depends_on_subdirs(ElementsKernel)
#===============================================================================
elements_depends_on_subdirs(EleFits EleFitsValidation)
elements_depends_on_subdirs(ElementsKernel)

#===============================================================================
# Add the find_package macro (a pure CMake command) here to locate the
# libraries.
# Examples:
#          find_package(CppUnit)
#===============================================================================

#===============================================================================
# Declare the library dependencies here
# Example:
#         elements_add_library(ElementsExamples src/Lib/*.cpp
#                     INCLUDE_DIRS Boost ElementsKernel
#                     LINK_LIBRARIES Boost ElementsKernel
#                     PUBLIC_HEADERS ElementsExamples)
#===============================================================================
elements_add_library(PhiBox src/lib/*.cpp
                     INCLUDE_DIRS EleFits ElementsKernel
                     LINK_LIBRARIES EleFits ElementsKernel
                     PUBLIC_HEADERS PhiBox)

#===============================================================================
# Declare the executables here
# Example:
# elements_add_executable(ElementsProgramExample src/Program/ProgramExample.cpp
#                        INCLUDE_DIRS Boost ElementsExamples
#                        LINK_LIBRARIES Boost ElementsExamples)
#===============================================================================
elements_add_executable(PhiTimeIntegration src/program/PhiTimeIntegration.cpp
                     INCLUDE_DIRS PhiBox
                     LINK_LIBRARIES PhiBox)

#===============================================================================
# Declare the Boost tests here
# Example:
# elements_add_unit_test(BoostClassExample tests/src/Boost/ClassExample_test.cpp
#                       EXECUTABLE BoostClassExample_test
#                       INCLUDE_DIRS ElementsExamples
#                       LINK_LIBRARIES ElementsExamples TYPE Boost)
#===============================================================================
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

#===============================================================================
# Use the following macro for python modules, scripts and aux files:
#  elements_install_python_modules()
#  elements_install_scripts()
#===============================================================================

#===============================================================================
# Declare the Python programs here
# Examples :
# elements_add_python_program(PythonProgramExample
#                             ElementsExamples.PythonProgramExample)
#===============================================================================

#===============================================================================
# Add the elements_install_conf_files macro
# Examples:
#          elements_install_conf_files()
#===============================================================================
elements_install_conf_files()
