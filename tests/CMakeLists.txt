ENABLE_TESTING()

set(Boost_USE_STATIC_LIBS ON)
set(BOOST_COMPONENTS unit_test_framework timer system)
SEARCH_FOR_BOOST()

include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories(${Boost_INCLUDE_DIRS})

macro(addUnitTest name)
  if(${UNIT_TESTS})
    add_executable(${name} ${name}.cpp ${HEADERS})
    target_link_libraries(${name} ${Boost_LIBRARIES})
    add_test(${name}Unit ${name})
  endif()
endmacro(addUnitTest)

macro(addBenchmark name)
  if(${BENCHMARKS})
    add_executable(${name} ${name}.cpp ${HEADERS})
    target_link_libraries(${name} ${Boost_LIBRARIES})
  endif()
endmacro(addBenchmark)

addUnitTest("VectorTest")
addUnitTest("InertiaTest")
addUnitTest("PTransformTest")
addUnitTest("AutoDiffTest")

addBenchmark("PTransformBench")

configure_file(sva_pickle_test.py ${CMAKE_CURRENT_BINARY_DIR}/sva_pickle_test.py)
add_test(svaPickleTestUnit sva_pickle_test.py)
