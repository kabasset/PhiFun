CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

elements_subdir(PhiFourier)

elements_depends_on_subdirs(LitlCore)

find_package(FFTW)

elements_add_library(PhiFourier src/lib/*.cpp
                     INCLUDE_DIRS LitlCore FFTW
                     LINK_LIBRARIES LitlCore FFTW
                     PUBLIC_HEADERS PhiFourier)

elements_add_unit_test(Dft tests/src/Dft_test.cpp 
                     EXECUTABLE PhiFourier_Dft_test
                     LINK_LIBRARIES PhiFourier
                     TYPE Boost)
elements_add_unit_test(DftBuffer tests/src/DftBuffer_test.cpp 
                     EXECUTABLE PhiFourier_DftBuffer_test
                     LINK_LIBRARIES PhiFourier
                     TYPE Boost)
elements_add_unit_test(DftMemory tests/src/DftMemory_test.cpp 
                     EXECUTABLE PhiFourier_DftMemory_test
                     LINK_LIBRARIES PhiFourier
                     TYPE Boost)
elements_add_unit_test(DftPlan tests/src/DftPlan_test.cpp 
                     EXECUTABLE PhiFourier_DftPlan_test
                     LINK_LIBRARIES PhiFourier
                     TYPE Boost)
