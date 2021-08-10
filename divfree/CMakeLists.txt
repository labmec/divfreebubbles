target_include_directories(divfree PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR})

set(MATERIAL_HEADERS
  TPZMatDivFreeBubbles.h
  TPZL2ProjectionCS.h
  TPZCompElKernelHdivBC.h)

set(MATERIAL_SRCS
  TPZMatDivFreeBubbles.cpp
  TPZL2ProjectionCS.cpp
  TPZCompElKernelHdivBC.cpp)


set(MESH_HEADERS
  TPZCompElKernelHdiv.h)

set(MESH_SRCS
  TPZCompElKernelHdiv.cpp)

target_sources(divfree PRIVATE
  ${MATERIAL_HEADERS} ${MATERIAL_SRCS}
  ${MESH_HEADERS} ${MESH_SRCS})
