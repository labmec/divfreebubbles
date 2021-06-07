cmake_minimum_required(VERSION 3.8)

project (DivFreeBubbles LANGUAGES C CXX)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# this option is useful for static analyzers
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Mac OS X specific code
IF(APPLE)
   add_definitions(-DMACOSX)
   set(MACOSX MACOSX)
   FIND_LIBRARY(ACCELERATE_LIB Accelerate)
   link_libraries(${ACCELERATE_LIB})
#   SET(OperatingSystem "Mac OS X")
ENDIF(APPLE)

# include_directories("${CMAKE_CURRENT_SOURCE_DIR}/headers")

# Finding the pz package
find_package(NeoPZ REQUIRED HINTS ../NeoPZ_install ../neopz_install)

# Adding pz files from this project
add_subdirectory(headers)
set(outsidepzset ${MATERIAL_FILES} ${MESH_FILES})

# creating a library with them and linking
add_library(pzfiles ${outsidepzset})
target_link_libraries(pzfiles PRIVATE NeoPZ::pz)
target_include_directories(pzfiles PRIVATE ${PZ_INCLUDE_DIRS})

# creating the project
add_subdirectory(src)

#organizing pz files of this project in IDE
source_group(Material FILES ${MATERIAL_FILES})
source_group(Mesh FILES ${MESH_FILES})
