# @file neopz/UnitTest_PZ/CMakeList.txt  -- First Cmake file for Unit Testing modules

enable_testing()
#in this file (and ONLY in this file) we will define CATCH_CONFIG_MAIN
add_library(test_library STATIC test_main.cpp)
if (APPLE)
  target_compile_definitions(test_library PUBLIC MACOSX)
endif()

target_link_libraries(test_library PUBLIC divfree Catch2::Catch2)

add_unit_test(edge_filter_test edge_filter_test.cpp)