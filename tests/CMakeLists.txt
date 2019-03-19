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

find_package(Boost REQUIRED COMPONENTS unit_test_framework timer)

macro(addUnitTest name)
  if(NOT ${DISABLE_TESTS})
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PUBLIC SpaceVecAlg Boost::unit_test_framework)
    add_test(${name}Unit ${name})
    # Adding a project configuration file (for MSVC only)
    GENERATE_MSVC_DOT_USER_FILE(${name})
  endif()
endmacro(addUnitTest)

macro(addBenchmark name)
  if(${BENCHMARKS})
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PUBLIC SpaceVecAlg Boost::unit_test_framework Boost::timer)
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
