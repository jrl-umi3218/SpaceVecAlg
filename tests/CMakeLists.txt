# Copyright 2012-2017 CNRS-UM LIRMM, CNRS-AIST JRL
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

ENABLE_TESTING()

set(BOOST_COMPONENTS unit_test_framework timer system)
SEARCH_FOR_BOOST()

include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

macro(addUnitTest name)
  if(NOT ${DISABLE_TESTS})
    add_executable(${name} ${name}.cpp ${HEADERS})
    if(NOT MSVC)
      # for msvc, we let boost auto-link do the work
      target_link_libraries(${name} ${Boost_LIBRARIES})
    endif(NOT MSVC)
    add_test(${name}Unit ${name})
    # Adding a project configuration file (for MSVC only)
    GENERATE_MSVC_DOT_USER_FILE(${name})
  endif()
endmacro(addUnitTest)

macro(addBenchmark name)
  if(${BENCHMARKS})
    add_executable(${name} ${name}.cpp ${HEADERS})
    if(NOT MSVC)
      # for msvc, we let boost auto-link do the work
      target_link_libraries(${name} ${Boost_LIBRARIES})
    endif(NOT MSVC)
    # Adding a project configuration file (for MSVC only)
    GENERATE_MSVC_DOT_USER_FILE(${name})
  endif()
endmacro(addBenchmark)

addUnitTest("VectorTest")
addUnitTest("InertiaTest")
addUnitTest("PTransformTest")
addUnitTest("AutoDiffTest")
addUnitTest("ConversionsTest")

addBenchmark("PTransformBench")
