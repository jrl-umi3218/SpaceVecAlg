#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

cmake_minimum_required(VERSION 3.1.3)

set(CMAKE_CXX_STANDARD 11)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(PROJECT_NAME SpaceVecAlg)
set(PROJECT_DESCRIPTION "Implementation of spatial vector algebra with the Eigen3 linear algebra library.")
set(PROJECT_URL "https://github.com/jrl-umi3218/SpaceVecAlg")
set(PROJECT_VERSION 1.1.0)
set(PROJECT_USE_CMAKE_EXPORT TRUE)

include(cmake/base.cmake)
include(cmake/cython/cython.cmake)
include(cmake/msvc-specific.cmake)

# Disable -Werror on Unix for now.
set(CXX_DISABLE_WERROR True)
project(${PROJECT_NAME} CXX)

option(BENCHMARKS "Generate benchmark." OFF)

add_project_dependency(Eigen3 REQUIRED)

# For MSVC, set local environment variable to enable finding the built dll
# of the main library when launching ctest with RUN_TESTS
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
