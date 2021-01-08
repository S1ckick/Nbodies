cmake_minimum_required(VERSION 3.17)
project(NBodies)

set(CMAKE_CXX_STANDARD 20 )
set(CMAKE_CXX_COMPILER /usr/local/bin/g++-10)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lquadmath -fext-numeric-literals -Wall -fexceptions -std=gnu++17 -g  -fpermissive ") #-mlong-double-128


include_directories("qd/2.3.22/include" "/usr/local/Cellar/boost/1.75.0_1/include")
link_directories("qd/2.3.22/lib" "/usr/local/Cellar/boost/1.75.0_1/lib")

add_executable(NBodies main.cpp Utils/vec.h Integration/methods.h Nbodies/nbodies.h Writer/writer.h Nbodies/energy.h Nbodies/summation.h Utils/helper.h)

target_link_libraries(NBodies qd)
target_link_libraries(NBodies qdmod)
target_link_libraries(NBodies qd_f_main)