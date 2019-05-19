cmake_minimum_required(VERSION 3.10.2)
project(ac370756)
set(CMAKE_CXX_STANDARD 14)
SET(CMAKE_CXX_COMPILER mpic++)

SET(CMAKE_CXX_FLAGS "-Wall -Wextra -g")

include_directories(.)

add_executable(matrixmul main.cpp
        utils.h utils.cpp
        parseInput.h parseInput.cpp
        CSCMatrix.h CSCMatrix.cpp
        const.h
        DenseMatrix.h DenseMatrix.cpp
        densematgen.h densematgen.cpp)
#target_link_libraries(matrixmul ${MPI_C_LIBRARIES})
#target_link_libraries(matrixmul libboost_program_options.a libboost_log.a)

FIND_PACKAGE( Boost 1.40 COMPONENTS program_options REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

TARGET_LINK_LIBRARIES( matrixmul LINK_PUBLIC ${Boost_LIBRARIES} )

configure_file(run.sh . COPYONLY)

add_custom_command(TARGET matrixmul
        POST_BUILD
        COMMAND /bin/sh /home/albert/HPC/MPI/ac370756/copy_tests.sh
        )
