cmake_minimum_required(VERSION 3.17)
project(NBodies)

set(CMAKE_CXX_STANDARD 20)

include_directories("/usr/local/Cellar/qd/2.3.22/include")


add_executable(NBodies main.cpp Utils/vec.h Integration/methods.h Nbodies/nbodies.h Writer/writer.h Nbodies/energy.h Nbodies/summation.h)

target_link_libraries(NBodies "/usr/local/Cellar/qd/2.3.22/lib/libqd.a")