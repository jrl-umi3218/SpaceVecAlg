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
# along with SpaceVecAlg.  If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 2.8)

INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/cpack.cmake)
INCLUDE(cmake/eigen.cmake)

SET(PROJECT_NAME SpaceVecAlg)
SET(PROJECT_DESCRIPTION "Implementation of spatial vector algebra with the Eigen3 linear algebra library.")
SET(PROJECT_URL "https://github.com/jrl-umi3218/SpaceVecAlg")

# Disable -Werror on Unix for now.
SET(CXX_DISABLE_WERROR True)
IF(NOT WIN32)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++0x -pedantic")
ENDIF()
SETUP_PROJECT()

option(PYTHON_BINDING "Generate python binding." OFF)
option(UNIT_TESTS "Generate unit tests." ON)
option(BENCHMARKS "Generate benchmark." OFF)

SET(Eigen_REQUIRED "eigen3 >= 3.2.0")
SEARCH_FOR_EIGEN()

add_subdirectory(src)

if(${UNIT_TESTS} OR ${BENCHMARKS})
  add_subdirectory(tests)
endif()

if(${PYTHON_BINDING})
  add_subdirectory(binding/python)
endif()

SETUP_PROJECT_FINALIZE()
SETUP_PROJECT_CPACK()

