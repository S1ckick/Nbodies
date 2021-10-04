cmake_minimum_required(VERSION 3.17)
project(NBodies)

set(CMAKE_CXX_STANDARD 20 )

include_directories("qd/2.3.22/include")
link_directories("qd/2.3.22/lib")

add_executable(NBodies main.cpp Utils/vec.h Integration/methods.h Nbodies/nbodies.h Writer/writer.h Nbodies/energy.h Nbodies/summation.h Utils/helper.h Integration/pointmasses.h)

target_link_libraries(NBodies qd)
target_link_libraries(NBodies qdmod)
target_link_libraries(NBodies qd_f_main)