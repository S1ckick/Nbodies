cmake_minimum_required(VERSION 3.17)
project(NBodies)

set(CMAKE_CXX_STANDARD 20)

add_executable(NBodies main.cpp Utils/vec.h Integration/methods.h Nbodies/nbodies.h Writer/writer.h Nbodies/energy.h Nbodies/summation.h)