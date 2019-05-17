cmake_minimum_required(VERSION 3.10.2)
project(ac370756)
set(CMAKE_CXX_STANDARD 14)
SET(CMAKE_CXX_COMPILER mpic++)

SET(CMAKE_CXX_FLAGS -lboost_program_options)

include_directories(.)

add_executable(matrixmul main.cpp)
target_link_libraries(matrixmul ${MPI_C_LIBRARIES})
target_link_libraries(matrixmul libboost_program_options.a)