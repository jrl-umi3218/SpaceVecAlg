# Copyright 2012-2016 CNRS-UM LIRMM, CNRS-AIST JRL
#
# This file is part of SpaceVecAlg.
#
# SpaceVecAlg is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SpaceVecAlg is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with SpaceVecAlg.  If not, see <http://www.gnu.or/licenses/>.

cmake_minimum_required(VERSION 3.1.3)

SET(PROJECT_NAME SpaceVecAlg)
SET(PROJECT_DESCRIPTION "Implementation of spatial vector algebra with the Eigen3 linear algebra library.")
SET(PROJECT_URL "https://github.com/jrl-umi3218/SpaceVecAlg")

INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/cpack.cmake)
INCLUDE(cmake/eigen.cmake)
INCLUDE(cmake/msvc-specific.cmake)

# Disable -Werror on Unix for now.
SET(CXX_DISABLE_WERROR True)
SETUP_PROJECT()
SET(PROJECT_VERSION 1.1.0)

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

SET(Eigen_REQUIRED "eigen3 >= 3.2.0")
SEARCH_FOR_EIGEN()
IF(NOT DEFINED Eigen_INCLUDE_DIRS)
  SET(Eigen_INCLUDE_DIRS "${_Eigen_INCLUDE_DIRS}")
ENDIF()

# For MSVC, set local environment variable to enable finding the built dll
# of the main library when launching ctest with RUN_TESTS
IF(MSVC)
  SET(CMAKE_MSVCIDE_RUN_PATH "\$(SolutionDir)/src/\$(Configuration)")
ENDIF(MSVC)

add_subdirectory(src)

if(NOT ${DISABLE_TESTS} OR ${BENCHMARKS})
  add_subdirectory(tests)
endif()

if(${PYTHON_BINDING})
  add_subdirectory(binding/python)
endif()

SETUP_PROJECT_FINALIZE()
SETUP_PROJECT_CPACK()
SETUP_PROJECT_PACKAGE_FINALIZE()
