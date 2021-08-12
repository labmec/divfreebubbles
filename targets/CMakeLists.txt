add_executable(dFreeBubbles main.cpp)
target_link_libraries(dFreeBubbles PRIVATE divfree)

add_executable(dFreeBubbles1el main_1element.cpp)
target_link_libraries(dFreeBubbles1el PRIVATE divfree)

add_executable(2frac main_2fractures.cpp)
target_link_libraries(2frac PRIVATE divfree)

add_executable(tetra_test tetra_test.cpp)
target_link_libraries(tetra_test PRIVATE divfree)

add_executable(semiHybrid semiHybrid.cpp)
target_link_libraries(semiHybrid PRIVATE divfree)