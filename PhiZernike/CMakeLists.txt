CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

elements_subdir(PhiZernike)

elements_depends_on_subdirs(LitlCore)
elements_depends_on_subdirs(ElementsKernel)

elements_add_library(PhiZernike src/lib/*.cpp
                     INCLUDE_DIRS LitlCore ElementsKernel
                     LINK_LIBRARIES LitlCore ElementsKernel
                     PUBLIC_HEADERS PhiZernike)

elements_add_unit_test(Chong tests/src/Chong_test.cpp 
                     EXECUTABLE PhiZernike_Chong_test
                     LINK_LIBRARIES PhiZernike
                     TYPE Boost)
elements_add_unit_test(Indexing tests/src/Indexing_test.cpp 
                     EXECUTABLE PhiZernike_Indexing_test
                     LINK_LIBRARIES PhiZernike
                     TYPE Boost)
elements_add_unit_test(Zernike tests/src/Zernike_test.cpp 
                     EXECUTABLE PhiZernike_Zernike_test
                     LINK_LIBRARIES PhiZernike
                     TYPE Boost)
