#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

cmake_minimum_required(VERSION 3.22)

find_package(jrl-cmakemodules 0.2.0 CONFIG QUIET)
if(jrl-cmakemodules_FOUND)
  message(
    STATUS "Using jrl-cmakemodules from CMake package: ${jrl-cmakemodules_DIR}")
else()
  message(STATUS "Using jrl-cmakemodules from source directory: cmake-v2")
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/cmake-v2 jrl-cmakemodules)
endif()

# Configure default binary and install directories, build type, etc.
jrl_configure_defaults()

project(
  SpaceVecAlg
  VERSION 1.2.7
  DESCRIPTION
    "Implementation of spatial vector algebra with the Eigen3 linear algebra library."
  HOMEPAGE_URL "https://github.com/jrl-umi3218/SpaceVecAlg"
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(PROJECT_USE_CMAKE_EXPORT TRUE)
set(PYTHON_BINDING_USE_PYTEST ON)

# Disable -Werror on Unix for now.
set(CXX_DISABLE_WERROR True)

# -------
# OPTIONS
# -------
jrl_option(BUILD_TESTING "Build testing." ON)
jrl_option(BENCHMARKS "Generate benchmark." OFF)
jrl_option(NANOBIND_BINDINGS "Build nanobind Python bindings" ON)
jrl_option(PYTHON_BINDING "Build nanobind Python bindings" OFF)

# ------------
# DEPENDENCIES
# ------------
jrl_find_package(Eigen3 REQUIRED)

# For MSVC, set local environment variable to enable finding the built dll of
# the main library when launching ctest with RUN_TESTS if(MSVC)
# set(CMAKE_MSVCIDE_RUN_PATH "\$(SolutionDir)/src/\$(Configuration)")
# endif(MSVC)

# FIXME(jrl-cmakemodules v2): when moving the code of this section to
# src/CMakeLists.txt it fails because it cannot find SpaceVecAlg in
# BUILDSYSTEM_TARGETS as BUILDSYSTEM_TARGETS is local to the directory being
# processed

# add_subdirectory(src)

set(HEADERS
    src/SpaceVecAlg/fwd.h
    src/SpaceVecAlg/MotionVec.h
    src/SpaceVecAlg/ForceVec.h
    src/SpaceVecAlg/ImpedanceVec.h
    src/SpaceVecAlg/AdmittanceVec.h
    src/SpaceVecAlg/PTransform.h
    src/SpaceVecAlg/RBInertia.h
    src/SpaceVecAlg/ABInertia.h
    src/SpaceVecAlg/EigenTypedef.h
    src/SpaceVecAlg/EigenUtility.h
    src/SpaceVecAlg/Operators.h
    src/SpaceVecAlg/MathFunc.h
    src/SpaceVecAlg/Conversions.h
    src/SpaceVecAlg/SpaceVecAlg)

add_library(SpaceVecAlg INTERFACE)
jrl_target_enforce_msvc_conformance(SpaceVecAlg PRIVATE)
target_link_libraries(SpaceVecAlg INTERFACE Eigen3::Eigen)
target_compile_features(SpaceVecAlg INTERFACE cxx_auto_type cxx_constexpr)
if(MSVC)
  target_compile_options(SpaceVecAlg INTERFACE /source-charset:utf-8)
endif()
target_include_directories(
  SpaceVecAlg INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
                        $<INSTALL_INTERFACE:include>)

jrl_target_headers(SpaceVecAlg PUBLIC HEADERS ${HEADERS})
set_target_properties(SpaceVecAlg PROPERTIES VERSION ${PROJECT_VERSION})

# Declare export components
jrl_add_export_component(NAME SpaceVecAlg TARGETS SpaceVecAlg)
# END FIXME

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

if(PYTHON_BINDING)
  add_subdirectory(binding/python)
endif()

# Export the package files (headers, targets, config files, etc.)
jrl_export_package()

# Print summaries
jrl_print_dependencies_summary()
jrl_print_options_summary()
