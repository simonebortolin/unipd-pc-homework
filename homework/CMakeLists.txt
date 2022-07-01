cmake_minimum_required(VERSION 3.16)
project(homework)

set(CMAKE_CXX_STANDARD 14)

add_executable(sequential main.cpp matrix.cpp matrix.h sequential.cpp sequential.h)

add_executable(parallel parallel.cpp matrix.cpp matrix.h sequential.cpp sequential.h)

add_executable(parallel2 parallel2.cpp matrix.cpp matrix.h sequential.cpp sequential.h)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
