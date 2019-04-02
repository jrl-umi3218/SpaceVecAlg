#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

enable_testing()

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
