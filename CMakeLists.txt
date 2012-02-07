cmake_minimum_required(VERSION 2.8)
project(SpaceVecAlg)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -std=c++0x -pedantic")

add_subdirectory(src)
add_subdirectory(tests)
