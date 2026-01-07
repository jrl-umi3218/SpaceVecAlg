#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

cmake_minimum_required(VERSION 3.22)

find_package(jrl-cmakemodules 0.2.0 CONFIG QUIET)
set(JRL_CMAKEMODULES_USE_V2
    ON
    CACHE BOOL "Use jrl-cmakemodules v2")
if(jrl-cmakemodules_FOUND)
  message(
    STATUS "Using jrl-cmakemodules from CMake package: ${jrl-cmakemodules_DIR}")
else()
  message(STATUS "Using jrl-cmakemodules from source directory: cmake-v2")
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/cmake-v2 jrl-cmakemodules)
endif()

project(
  SpaceVecAlg
  VERSION 1.2.7
  DESCRIPTION
    "Implementation of spatial vector algebra with the Eigen3 linear algebra library."
  HOMEPAGE_URL "https://github.com/jrl-umi3218/SpaceVecAlg"
  LANGUAGES CXX)

# Configure default binary and install directories, build type, etc.
jrl_configure_defaults()

set(CMAKE_CXX_STANDARD 17)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# -------
# OPTIONS
# -------
jrl_option(BUILD_TESTING "Build testing." ON)
jrl_option(BENCHMARKS "Generate benchmark." OFF)
jrl_option(NANOBIND_BINDINGS "Build nanobind Python bindings" ON)
# Disable old python bindings jrl_option(PYTHON_BINDING "Build cython Python
# bindings" OFF)

# ------------
# DEPENDENCIES
# ------------
jrl_find_package(Eigen3 REQUIRED)

# For MSVC, set local environment variable to enable finding the built dll of
# the main library when launching ctest with RUN_TESTS
if(MSVC)
  set(CMAKE_MSVCIDE_RUN_PATH "\$(SolutionDir)/src/\$(Configuration)")
endif(MSVC)

add_subdirectory(src)

if(BUILD_TESTING OR BENCHMARKS)
  enable_testing()
  add_subdirectory(tests)
endif()

if(NANOBIND_BINDINGS)
  message(STATUS "[${PROJECT_NAME}] Building Python bindings")
  jrl_find_python(3.8 REQUIRED COMPONENTS Interpreter Development.Module)
  jrl_find_nanobind(2.4.0 CONFIG REQUIRED)
  add_subdirectory(binding/nanobind)
endif()

# Disable cython bindings if(PYTHON_BINDING) add_subdirectory(binding/python)
# endif()

# Export the package files (headers, targets, config files, etc.)
jrl_export_package()

# Print summaries
jrl_print_dependencies_summary()
jrl_print_options_summary()
