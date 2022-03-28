CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

#===============================================================================
# Load macros and functions for Elements-based projects
#===============================================================================

find_package(ElementsProject)

set(ELEMENTS_PARALLEL ON
    CACHE STRING "Enable parallelization support with OpenMP"
    FORCE)

#===============================================================================
# Declare project name and version
# Example with dependency:
#                         elements_project(MyProject 1.0 USE Element 3.9)
#===============================================================================

elements_project(PhiFun 0.9 USE EleFits 5.1 DESCRIPTION "Fun with Phi!")
