# Set the minimum CMake version and policies for highest tested version
cmake_minimum_required(VERSION 3.15...4.0)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(SpaceVecAlg)

# Warn if the user invokes CMake directly
if(NOT SKBUILD)
  message(
    WARNING
      "\
  This CMake file is meant to be executed using 'scikit-build-core'.
  Running it directly will almost certainly not produce the desired
  result. If you are a user trying to install this package, use the
  command below, which will install all necessary build dependencies,
  compile the package in an isolated environment, and then install it.
  =====================================================================
   $ pip install .
  =====================================================================
  If you are a software developer, and this is your own package, then
  it is usually much more efficient to install the build dependencies
  in your environment once and use the following command that avoids
  a costly creation of a new virtual environment at every compilation:
  =====================================================================
   $ pip install nanobind scikit-build-core[pyproject]
   $ pip install --no-build-isolation -ve .
  =====================================================================
  You may optionally add -Ceditable.rebuild=true to auto-rebuild when
  the package is imported. Otherwise, you need to rerun the above
  after editing C++ files.")
endif()

# Try to import all Python components potentially needed by nanobind
find_package(
  Python 3.8 REQUIRED
  COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule)

# Import nanobind through CMake's find_package mechanism
find_package(nanobind CONFIG REQUIRED)

# We are now ready to compile the actual extension module
nanobind_add_module(
  # Name of the extension
  _sva
  # Target the stable ABI for Python 3.12+, which reduces the number of binary
  # wheels that must be built. This does nothing on older Python versions
  STABLE_ABI
  # Build libnanobind statically and merge it into the extension (which itself
  # remains a shared library)
  #
  # If your project builds multiple extensions, you can replace this flag by
  # NB_SHARED to conserve space by reusing a shared libnanobind across libraries
  NB_STATIC
  # Source code goes here
  src/sva_module/sva_module.cpp
  src/sva_module/bind_PTransformd.cpp)

find_package(SpaceVecAlg REQUIRED)
target_link_libraries(_sva PUBLIC SpaceVecAlg::SpaceVecAlg)
target_include_directories(
  _sva PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
# Install directive for scikit-build-core
install(TARGETS _sva LIBRARY DESTINATION sva)

# Display output of stub generation
set(CMAKE_EXECUTE_PROCESS_COMMAND_ECHO STDOUT)

nanobind_add_stub(
  _sva_stub
  MODULE
  sva
  OUTPUT
  ${SKBUILD_PLATLIB_DIR}/sva/__init__.pyi
  MARKER_FILE
  ${SKBUILD_PLATLIB_DIR}/sva/py.typed
  PYTHON_PATH
  ${SKBUILD_PLATLIB_DIR}/sva/
  DEPENDS
  _sva)
