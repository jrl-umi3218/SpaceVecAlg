#
# Copyright 2012-2025 CNRS-UM LIRMM, CNRS-AIST JRL
#

jrl_find_package(Boost CONFIG REQUIRED COMPONENTS unit_test_framework timer)

macro(addUnitTest name)
  if(${BUILD_TESTING})
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PUBLIC SpaceVecAlg Boost::unit_test_framework)
    set_target_properties(${name} PROPERTIES FOLDER "tests")
    add_test(
      NAME ${name}Unit
      COMMAND ${name}
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
    # Adding a project configuration file (for MSVC only)
    # generate_msvc_dot_user_file(${name})
  endif()
endmacro(addUnitTest)

macro(addBenchmark name)
  if(${BENCHMARKS})
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PUBLIC SpaceVecAlg Boost::unit_test_framework
                                         Boost::timer)
    # Adding a project configuration file (for MSVC only)
    # generate_msvc_dot_user_file(${name})
  endif()
endmacro(addBenchmark)

addunittest("VectorTest")
addunittest("InertiaTest")
addunittest("PTransformTest")
addunittest("AutoDiffTest")
addunittest("ConversionsTest")
addunittest("LogDiffTest")

addbenchmark("PTransformBench")
