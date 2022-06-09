target_include_directories(divfree PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR})

set(MATERIAL_HEADERS
  TPZMatDivFreeBubbles.h
  TPZMatCurlDotCurl.h
  TPZMixedDarcyFlowHybrid.h
  )

set(MATERIAL_SRCS
  TPZMatDivFreeBubbles.cpp
  TPZMatCurlDotCurl.cpp
  TPZMixedDarcyFlowHybrid.cpp
  )


set(MESH_HEADERS
    TPZCompElHDivDuplConnects.h
    TPZCompElHDivDuplConnectsBound.h
  )

set(MESH_SRCS
    TPZCompElHDivDuplConnects.cpp
    TPZCompElHDivDuplConnectsBound.cpp
  )

# let us create a config file with the directory containing the meshes
# create folder with meshes
set(MESHESDIR ${PROJECT_SOURCE_DIR}/mesh)
set (CONFIG_H "${CMAKE_CURRENT_SOURCE_DIR}/divfree_config.h")
configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/divfree_config.h.in" ${CONFIG_H})

set(UTIL_HEADERS
  divfree_config.h
  TPZHDivApproxSpaceCreator.h
  TPZKernelHdivHybridizer.h
  TPZKernelHdivUtils.h
  TPZSparseMatRed.h
  TPZMatRedSolver.h
  TPZCompElConstFluxHybrid.h
  )

set(UTIL_SRCS
  TPZKernelHdivHybridizer.cpp
  TPZHDivApproxSpaceCreator.cpp
  TPZKernelHdivUtils.cpp
  TPZSparseMatRed.cpp
  TPZMatRedSolver.cpp
  TPZCompElConstFluxHybrid.cpp
  )

target_sources(divfree PRIVATE
  ${MATERIAL_HEADERS} ${MATERIAL_SRCS}
  ${MESH_HEADERS} ${MESH_SRCS}
  ${UTIL_HEADERS} ${UTIL_SRCS})