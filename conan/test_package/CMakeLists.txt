cmake_minimum_required(VERSION 3.1.3)
set(CMAKE_CXX_STANDARD 17)
project(PackageTest CXX)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

add_executable(example example.cpp)
target_link_libraries(example ${CONAN_LIBS})
