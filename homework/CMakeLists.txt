cmake_minimum_required(VERSION 3.16)
project(homework)

set(CMAKE_CXX_STANDARD 14)

add_executable(homework main.cpp matrix.cpp matrix.h sequential.cpp sequential.h)

add_executable(parallel parallel.cpp matrix.cpp matrix.h)