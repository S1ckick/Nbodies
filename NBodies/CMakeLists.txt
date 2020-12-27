cmake_minimum_required(VERSION 3.17)
project(NBodies)

set(CMAKE_CXX_STANDARD 20)

add_executable(NBodies main.cpp Math/vec.h Integration/integration.h Integration/RungeKutta4.h Bodies/bodies.h)