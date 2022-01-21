cmake_minimum_required(VERSION 3.14)
project(NBodies LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

include_directories("qd/2.3.22/include")



add_executable(NBodies main.cpp)


target_link_libraries(NBodies qd qdmod qd_f_main)