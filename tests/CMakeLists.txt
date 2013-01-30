ENABLE_TESTING()

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost COMPONENTS unit_test_framework REQUIRED)

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

