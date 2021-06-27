# add_pz_target(
#     NAME dFreeBubbles
#     SOURCES main.cpp)

# add_pz_target(
#     NAME dFreeBubbles1el
#     SOURCES main_1element.cpp)

add_executable(dFreeBubbles main.cpp)
target_link_libraries(dFreeBubbles PRIVATE pzfiles NeoPZ::pz)

add_executable(dFreeBubbles1el main_1element.cpp)
target_link_libraries(dFreeBubbles1el PRIVATE pzfiles NeoPZ::pz)

add_executable(2frac main_2fractures.cpp)
target_link_libraries(2frac PRIVATE pzfiles NeoPZ::pz)