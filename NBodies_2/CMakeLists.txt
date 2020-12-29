cmake_minimum_required(VERSION 3.17)
project(NBodies_2)

set(CMAKE_CXX_STANDARD 20)

add_executable(NBodies_2 main.cpp Utils/vec.h Integration/methods.h Nbodies/nbodies.h Writer/writer.h Nbodies/integration.h)