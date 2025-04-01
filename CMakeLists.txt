cmake_minimum_required(VERSION 3.16)

project(libdgp VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 17)

list(PREPEND CMAKE_MODULE_PATH 
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake
)

option(LIBDGP_BUILD_TUTORIALS  "Build libdgp tutorial" ON)
option(LIBDGP_USE_STATIC_LIBRARY "Use libdgp as static library" OFF)

include(dgp_add_tutorial)

include(libigl)
include(libdgp)
#include(geometry-central)

# Enable the target of libigl
igl_include(glfw)
igl_include(imgui)

if(LIBDGP_BUILD_TUTORIALS)
    add_subdirectory(tutorial)
endif()
