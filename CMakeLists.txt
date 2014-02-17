cmake_minimum_required(VERSION 2.8)

INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/cpack.cmake)
INCLUDE(cmake/eigen.cmake)

SET(PROJECT_NAME SpaceVecAlg)
SET(PROJECT_DESCRIPTION "Implementation of spatial vector algebra with the Eigen3 linear algebra library.")
SET(PROJECT_URL "https://github.com/jorisv/SpaceVecAlg")

# Disable -Werror on Unix for now.
SET(CXX_DISABLE_WERROR True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++0x -pedantic")
SETUP_PROJECT()

option(PYTHON_BINDING "Generate python binding." ON)
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

