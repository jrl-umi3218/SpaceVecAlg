ENABLE_TESTING()

set(Boost_USE_STATIC_LIBS ON)
set(BOOST_COMPONENTS unit_test_framework)
SEARCH_FOR_BOOST()

include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories(${Boost_INCLUDE_DIRS})

add_executable(VectorTest VectorTest.cpp)
target_link_libraries(VectorTest ${Boost_LIBRARIES})
add_test(VectorTestUnit VectorTest)

add_executable(InertiaTest InertiaTest.cpp)
target_link_libraries(InertiaTest ${Boost_LIBRARIES})
add_test(InertiaTestUnit InertiaTest)

add_executable(PTransformTest PTransformTest.cpp)
target_link_libraries(PTransformTest ${Boost_LIBRARIES})
add_test(PTransformTestUnit PTransformTest)

add_executable(AutoDiffTest AutoDiffTest.cpp)
target_link_libraries(AutoDiffTest ${Boost_LIBRARIES})
add_test(AutoDiffTestUnit AutoDiffTest)

