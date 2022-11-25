add_executable(dFreeBubbles main.cpp)
target_link_libraries(dFreeBubbles PRIVATE divfree)

add_executable(dFreeBubbles1el main_1element.cpp)
target_link_libraries(dFreeBubbles1el PRIVATE divfree)

add_executable(2frac main_2fractures.cpp)
target_link_libraries(2frac PRIVATE divfree)

add_executable(semiHybrid semiHybrid.cpp)
target_link_libraries(semiHybrid PRIVATE divfree)

add_executable(3d main3d.cpp)
target_link_libraries(3d PRIVATE divfree)

add_executable(SHIterative SHIterative.cpp)
target_link_libraries(SHIterative PRIVATE divfree)

add_executable(fiveSpot fiveSpot.cpp)
target_link_libraries(fiveSpot PRIVATE divfree)

add_executable(3d_hdiv main3d_HDIV.cpp)
target_link_libraries(3d_hdiv PRIVATE divfree)

add_executable(MHM_HDivConstant main_MHM_HDivConstant.cpp)
target_link_libraries(MHM_HDivConstant PRIVATE divfree)

add_executable(dupl_connects2 dupl_connects.cpp)
target_link_libraries(dupl_connects2 PRIVATE divfree)
