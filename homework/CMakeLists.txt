cmake_minimum_required(VERSION 3.11)
project(homework)

set(CMAKE_CXX_STANDARD 14)
set(GCC_COVERAGE_COMPILE_FLAGS "-Wall -pedantic -lm -O3 -funroll-loops")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

add_executable(sequential main.cpp matrix.cpp matrix.h sequential.cpp sequential.h basic_matrix.h)

add_executable(parallel parallel.cpp matrix.cpp matrix.h sequential.cpp sequential.h matrix_matrix.h basic_matrix.h)


# OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(sequential PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(parallel PUBLIC OpenMP::OpenMP_CXX)
endif()

# MPI
find_package(MPI REQUIRED)
target_link_libraries(sequential PUBLIC MPI::MPI_CXX)
target_link_libraries(parallel PUBLIC MPI::MPI_CXX)


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
