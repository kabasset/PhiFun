CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

find_package(ElementsProject)

# set(ELEMENTS_PARALLEL ON
#     CACHE STRING "Enable parallelization support with OpenMP"
#     FORCE)

set(CMAKE_CXX_FLAGS_RELEASE "-flto -O3"
    CACHE STRING "Enable -o3 for release builds."
    FORCE)
set(CMAKE_C_FLAGS_RELEASE "-O3"
    CACHE STRING "Enable -o3 for release builds."
    FORCE)

elements_project(PhiFun 0.9 USE Litl 1.0 EleFits 5.1 DESCRIPTION "Fun with Phi!")
