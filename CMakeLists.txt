cmake_minimum_required(VERSION 3.8)

project (DivFreeBubbles LANGUAGES CXX)
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


# Finding the pz package
find_package(NeoPZ REQUIRED)

# adding the project
add_subdirectory(src)