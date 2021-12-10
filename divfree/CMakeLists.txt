target_include_directories(divfree PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR})

set(MATERIAL_HEADERS
  TPZMatDivFreeBubbles.h
  TPZL2ProjectionCS.h
  TPZCompElKernelHdivBC.h
  TPZCompElKernelHdivBC3D.h
  TPZMixedDarcyFlowHybrid.h
  TPZMatCurlDotCurl.h)

set(MATERIAL_SRCS
  TPZMatDivFreeBubbles.cpp
  TPZL2ProjectionCS.cpp
  TPZCompElKernelHdivBC.cpp
  TPZCompElKernelHdivBC3D.cpp
  TPZMixedDarcyFlowHybrid.cpp
  TPZMatCurlDotCurl.cpp)


set(MESH_HEADERS
  TPZCompElKernelHdiv.h
  TPZCompElHDivConstant.h
  TPZCompElHDivConstantBC.h
  TPZCompElKernelHDiv3D.h
  TPZCompElHDivConstant3D.h
  TPZCompElHCurlNoGrads.h)

set(MESH_SRCS
  TPZCompElKernelHdiv.cpp
  TPZCompElHDivConstant.cpp
  TPZCompElHDivConstantBC.cpp
  TPZCompElKernelHDiv3D.cpp
  TPZCompElHDivConstant3D.cpp
  TPZCompElHCurlNoGrads.cpp)

# let us create a config file with the directory containing the meshes
# create folder with meshes
set(MESHESDIR ${PROJECT_SOURCE_DIR}/mesh)
set (CONFIG_H "${CMAKE_CURRENT_SOURCE_DIR}/divfree_config.h")
configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/divfree_config.h.in" ${CONFIG_H})

set(UTIL_HEADERS
  divfree_config.h
  TPZApproxSpaceKernelHdiv.h
  TPZKernelHdivHybridizer.h
  TPZKernelHdivUtils.h
  TPZHCurlEquationFilter.h)

set(UTIL_SRCS
  TPZKernelHdivHybridizer.cpp
  TPZApproxSpaceKernelHdiv.cpp
  TPZKernelHdivUtils.cpp
  TPZHCurlEquationFilter.cpp)

target_sources(divfree PRIVATE
  ${MATERIAL_HEADERS} ${MATERIAL_SRCS}
  ${MESH_HEADERS} ${MESH_SRCS}
  ${UTIL_HEADERS} ${UTIL_SRCS})