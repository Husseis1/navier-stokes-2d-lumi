cmake_minimum_required(VERSION 3.10)
project(NavierStokes CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(MPI REQUIRED)

include_directories(${MPI_CXX_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/part1)

add_executable(part1_ns main.cc ns2d.cc)
target_link_libraries(part1_ns PUBLIC MPI::MPI_CXX)