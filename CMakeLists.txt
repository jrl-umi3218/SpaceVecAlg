#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

cmake_minimum_required(VERSION 3.22)

set(CMAKE_CXX_STANDARD 17)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(PROJECT_NAME SpaceVecAlg)
set(PROJECT_DESCRIPTION
    "Implementation of spatial vector algebra with the Eigen3 linear algebra library."
)
set(PROJECT_URL "https://github.com/jrl-umi3218/SpaceVecAlg")
set(PROJECT_VERSION 1.2.8)
set(PROJECT_USE_CMAKE_EXPORT TRUE)
set(PYTHON_BINDING_USE_PYTEST ON)

find_package(jrl-cmakemodules CONFIG)
get_property(
JRL_CMAKE_MODULES
TARGET jrl-cmakemodules::jrl-cmakemodules
PROPERTY INTERFACE_INCLUDE_DIRECTORIES)

set(JRL_CMAKEMODULES_USE_V2 ON CACHE BOOL "Use jrl-cmakemodules v2")
find_package(jrl-cmakemodules 0.2.0 CONFIG QUIET)
if(jrl-cmakemodules_FOUND)
else()
  message(STATUS "Using jrl-cmakemodules from source directory: cmake-v2")
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/cmake-v2 jrl-cmakemodules)
endif()

include(${JRL_CMAKE_MODULES}/base.cmake)
include(${JRL_CMAKE_MODULES}/cython/cython.cmake)
include(${JRL_CMAKE_MODULES}/msvc-specific.cmake)

# Disable -Werror on Unix for now.
set(CXX_DISABLE_WERROR True)
project(${PROJECT_NAME} CXX)

option(BENCHMARKS "Generate benchmark." OFF)

find_package(Eigen3 QUIET NO_CMAKE_PACKAGE_REGISTRY)
if(Eigen3_FOUND)
  add_project_dependency(Eigen3 REQUIRED NO_CMAKE_PACKAGE_REGISTRY)
else()
  add_project_dependency(Eigen3 MODULE REQUIRED)
endif()

# For MSVC, set local environment variable to enable finding the built dll of
# the main library when launching ctest with RUN_TESTS
if(MSVC)
  set(CMAKE_MSVCIDE_RUN_PATH "\$(SolutionDir)/src/\$(Configuration)")
endif(MSVC)

add_subdirectory(src)

if(${BUILD_TESTING} OR ${BENCHMARKS})
  add_subdirectory(tests)
endif()

if(${PYTHON_BINDING})
  add_subdirectory(binding/python)
endif()
