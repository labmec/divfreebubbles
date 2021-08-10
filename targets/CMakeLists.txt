add_executable(dFreeBubbles main.cpp)
target_link_libraries(dFreeBubbles PRIVATE divfreelib)

add_executable(dFreeBubbles1el main_1element.cpp)
target_link_libraries(dFreeBubbles1el PRIVATE divfreelib)

add_executable(2frac main_2fractures.cpp)
target_link_libraries(2frac PRIVATE divfreelib)