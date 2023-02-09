# @file neopz/UnitTest_PZ/CMakeList.txt  -- First Cmake file for Unit Testing modules

enable_testing()
#in this file (and ONLY in this file) we will define CATCH_CONFIG_MAIN
add_library(test_library STATIC test_main.cpp)
if (APPLE)
  target_compile_definitions(test_library PUBLIC MACOSX)
endif()

target_link_libraries(test_library PUBLIC divfree Catch2::Catch2)

add_unit_test(edge_filter_test edge_filter_test.cpp)

add_unit_test(hybridization_test hybridization_test.cpp)

add_unit_test(dupl_connects dupl_connects.cpp)
add_unit_test(dMatRed doubleMatRed.cpp)

add_unit_test(mixed_elasticity mixed_elasticity.cpp)
add_unit_test(Yotov Yotov.cpp)

configure_file(InputData/spe_perm.dat InputData/spe_perm.dat COPYONLY)

add_unit_test(spe10 spe10.cpp)

add_unit_test(PGD PGD.cpp)
