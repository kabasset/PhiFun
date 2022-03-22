CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

#===============================================================================
# Load elements_subdir macro here
# Examples:
#   For declaring a project module:
#         elements_subdir(ElementsExamples)
#===============================================================================
elements_subdir(PhiDuffieux)

#===============================================================================
# Load elements_depends_on_subdirs macro here
#   For creating a dependency onto an other accessible module
#         elements_depends_on_subdirs(ElementsKernel)
#===============================================================================
elements_depends_on_subdirs(PhiZernike PhiFourier PhiBox)
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
elements_add_library(PhiDuffieux src/lib/*.cpp
                     INCLUDE_DIRS PhiZernike PhiFourier PhiBox ElementsKernel
                     LINK_LIBRARIES PhiZernike PhiFourier PhiBox ElementsKernel
                     PUBLIC_HEADERS PhiDuffieux)

#===============================================================================
# Declare the executables here
# Example:
# elements_add_executable(ElementsProgramExample src/Program/ProgramExample.cpp
#                        INCLUDE_DIRS Boost ElementsExamples
#                        LINK_LIBRARIES Boost ElementsExamples)
#===============================================================================
elements_add_executable(PhiEstimatePsf src/program/PhiEstimatePsf.cpp
                     INCLUDE_DIRS PhiDuffieux
                     LINK_LIBRARIES PhiDuffieux)
elements_add_executable(PhiGeneratePsf src/program/PhiGeneratePsf.cpp
                     INCLUDE_DIRS PhiDuffieux
                     LINK_LIBRARIES PhiDuffieux)

#===============================================================================
# Declare the Boost tests here
# Example:
# elements_add_unit_test(BoostClassExample tests/src/Boost/ClassExample_test.cpp
#                       EXECUTABLE BoostClassExample_test
#                       INCLUDE_DIRS ElementsExamples
#                       LINK_LIBRARIES ElementsExamples TYPE Boost)
#===============================================================================
elements_add_unit_test(MonochromaticOptics tests/src/MonochromaticOptics_test.cpp 
                     EXECUTABLE PhiDuffieux_MonochromaticOptics_test
                     LINK_LIBRARIES PhiDuffieux
                     TYPE Boost)
elements_add_unit_test(MonochromaticSystem tests/src/MonochromaticSystem_test.cpp 
                     EXECUTABLE PhiDuffieux_MonochromaticSystem_test
                     LINK_LIBRARIES PhiDuffieux
                     TYPE Boost)
elements_add_unit_test(BroadbandSystem tests/src/BroadbandSystem_test.cpp 
                     EXECUTABLE PhiDuffieux_BroadbandSystem_test
                     LINK_LIBRARIES PhiDuffieux
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
