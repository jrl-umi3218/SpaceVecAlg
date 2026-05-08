#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

cmake_minimum_required(VERSION 3.22)

find_package(jrl-cmakemodules 0.2.0 CONFIG QUIET)

project(
  SpaceVecAlg
  VERSION 1.2.7
  DESCRIPTION
    "Implementation of spatial vector algebra with the Eigen3 linear algebra library."
  HOMEPAGE_URL "https://github.com/jrl-umi3218/SpaceVecAlg"
  LANGUAGES CXX)

include(cmake/get-jrl-cmakemodules.cmake)

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
jrl_option(
  NANOBIND_DOCUMENTATION
  "Build nanobind documentation"
  OFF
  CONDITION
  NANOBIND_BINDINGS
  FALLBACK
  OFF)
jrl_option(BUILD_DOCUMENTATION "Whether doxygen documentation should be built"
           OFF)
jrl_option(
  INSTALL_DOCUMENTATION
  "Generate and install the doxygen documentation"
  OFF
  CONDITION
  BUILD_DOCUMENTATION
  FALLBACK
  OFF)

# Disable old python bindings jrl_option(PYTHON_BINDING "Build cython Python
# bindings" OFF)
jrl_print_options_summary()

# ------------
# DEPENDENCIES
# ------------
jrl_find_package(Eigen3 REQUIRED)
if(BUILD_DOCUMENTATION)
  jrl_find_package(Doxygen REQUIRED COMPONENTS dot)
endif()

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

if(BUILD_DOCUMENTATION)
  # Configure Doxygen
  set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/share/doc/sva/doxygen)
  set(DOXYGEN_USE_MDFILE_AS_MAINPAGE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
  doxygen_add_docs(doc ${PROJECT_SOURCE_DIR} ALL
                   COMMENT "Generating API documentation with Doxygen")

  if(INSTALL_DOCUMENTATION)
    install(
      DIRECTORY ${CMAKE_BINARY_DIR}/share/doc/sva/doxygen/html
      DESTINATION ${CMAKE_INSTALL_PREFIX}/share/doc/sva/doxygen
      COMPONENT doc)
  endif()
endif()

if(NANOBIND_BINDINGS)
  message(STATUS "[${PROJECT_NAME}] Building Python bindings")
  jrl_find_python(3.8 REQUIRED COMPONENTS Interpreter Development.Module)
  jrl_find_nanobind(2.4.0 CONFIG REQUIRED)
  jrl_find_package(nanoeigenpy CONFIG REQUIRED)
  if(NANOBIND_DOCUMENTATION)
    include(cmake/get-sphinx-cmake.cmake)
    add_subdirectory(binding/nanobind)
  endif()
endif()

# Disable cython bindings if(PYTHON_BINDING) add_subdirectory(binding/python)
# endif()

# Export the package files (headers, targets, config files, etc.)
jrl_export_package()

# Print summaries
jrl_print_dependencies_summary()
