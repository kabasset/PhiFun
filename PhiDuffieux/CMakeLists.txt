CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

elements_subdir(PhiDuffieux)

elements_depends_on_subdirs(PhiZernike PhiFourier PhiBox)
elements_depends_on_subdirs(EleFits)
elements_depends_on_subdirs(ElementsKernel)

# find_package(OpenMP)

elements_add_library(PhiDuffieux src/lib/*.cpp
                     INCLUDE_DIRS PhiZernike PhiFourier PhiBox ElementsKernel
                     LINK_LIBRARIES PhiZernike PhiFourier PhiBox ElementsKernel
                     PUBLIC_HEADERS PhiDuffieux)

elements_add_executable(PhiEstimatePsf src/program/PhiEstimatePsf.cpp
                     INCLUDE_DIRS PhiDuffieux EleFits # OpenMP
                     LINK_LIBRARIES PhiDuffieux EleFits) # OpenMP
elements_add_executable(PhiGeneratePsf src/program/PhiGeneratePsf.cpp
                     INCLUDE_DIRS PhiDuffieux EleFits
                     LINK_LIBRARIES PhiDuffieux EleFits)

elements_add_unit_test(BroadbandSystem tests/src/BroadbandSystem_test.cpp 
                     EXECUTABLE PhiDuffieux_BroadbandSystem_test
                     LINK_LIBRARIES PhiDuffieux
                     TYPE Boost)
elements_add_unit_test(MonochromaticOptics tests/src/MonochromaticOptics_test.cpp 
                     EXECUTABLE PhiDuffieux_MonochromaticOptics_test
                     LINK_LIBRARIES PhiDuffieux
                     TYPE Boost)
elements_add_unit_test(MonochromaticSystem tests/src/MonochromaticSystem_test.cpp 
                     EXECUTABLE PhiDuffieux_MonochromaticSystem_test
                     LINK_LIBRARIES PhiDuffieux
                     TYPE Boost)

elements_install_conf_files()
