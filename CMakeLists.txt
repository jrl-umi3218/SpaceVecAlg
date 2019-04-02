#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

cmake_minimum_required(VERSION 3.1.3)

set(PROJECT_NAME SpaceVecAlg)
set(PROJECT_DESCRIPTION "Implementation of spatial vector algebra with the Eigen3 linear algebra library.")
set(PROJECT_URL "https://github.com/jrl-umi3218/SpaceVecAlg")

include(cmake/base.cmake)
include(cmake/boost.cmake)
include(cmake/cpack.cmake)
include(cmake/msvc-specific.cmake)

# Disable -Werror on Unix for now.
set(CXX_DISABLE_WERROR True)
setup_project()
set(PROJECT_VERSION 1.1.0)

option(PYTHON_BINDING "Generate python binding." ON)
if(WIN32)
  set(PYTHON_BINDING_USER_INSTALL_DEFAULT ON)
else()
  set(PYTHON_BINDING_USER_INSTALL_DEFAULT OFF)
endif()
option(PYTHON_BINDING_USER_INSTALL "Install the Python bindings in user space" ${PYTHON_BINDING_USER_INSTALL_DEFAULT})
option(PYTHON_BINDING_FORCE_PYTHON3 "Use pip3/python3 instead of pip/python" OFF)
option(DISABLE_TESTS "Disable unit tests." OFF)
option(BENCHMARKS "Generate benchmark." OFF)

find_package(Eigen3 REQUIRED)

# For MSVC, set local environment variable to enable finding the built dll
# of the main library when launching ctest with RUN_TESTS
if(MSVC)
  set(CMAKE_MSVCIDE_RUN_PATH "\$(SolutionDir)/src/\$(Configuration)")
endif(MSVC)

add_subdirectory(src)

if(NOT ${DISABLE_TESTS} OR ${BENCHMARKS})
  add_subdirectory(tests)
endif()

if(${PYTHON_BINDING})
  add_subdirectory(binding/python)
endif()

setup_project_finalize()
setup_project_cpack()
setup_project_package_finalize()
