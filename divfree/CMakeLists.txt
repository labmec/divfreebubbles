target_include_directories(divfree PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR})

set(MATERIAL_HEADERS
  TPZMatDivFreeBubbles.h
  TPZL2ProjectionCS.h
  TPZCompElKernelHdivBC.h
  TPZMixedDarcyFlowHybrid.h
  TPZMatCurlDotCurl.h)

set(MATERIAL_SRCS
  TPZMatDivFreeBubbles.cpp
  TPZL2ProjectionCS.cpp
  TPZCompElKernelHdivBC.cpp
  TPZMixedDarcyFlowHybrid.cpp
  TPZMatCurlDotCurl.cpp)


set(MESH_HEADERS
  TPZCompElKernelHdiv.h)

set(MESH_SRCS
  TPZCompElKernelHdiv.cpp)

# let us create a config file with the directory containing the meshes
# create folder with meshes
set(MESHESDIR ${PROJECT_SOURCE_DIR}/mesh)
set (CONFIG_H "${CMAKE_CURRENT_SOURCE_DIR}/divfree_config.h")
configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/divfree_config.h.in" ${CONFIG_H})

set(UTIL_HEADERS
  divfree_config.h
  TPZKernelHdivHybridizer.h
  TPZKernelHdivUtils.h)

set(UTIL_SRCS
  TPZKernelHdivHybridizer.cpp
  TPZKernelHdivUtils.cpp)

target_sources(divfree PRIVATE
  ${MATERIAL_HEADERS} ${MATERIAL_SRCS}
  ${MESH_HEADERS} ${MESH_SRCS}
  ${UTIL_HEADERS} ${UTIL_SRCS})