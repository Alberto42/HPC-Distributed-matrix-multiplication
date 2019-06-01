cmake_minimum_required(VERSION 3.10.2)
project(ac370756)
set(CMAKE_CXX_STANDARD 14)
SET(CMAKE_CXX_COMPILER mpic++)

SET(CMAKE_CXX_FLAGS "-Wall -Wextra -O3 -g")

include_directories(.)


file(GLOB_RECURSE Foo_SOURCES "src/*.cpp")
file(GLOB_RECURSE Foo_HEADERS "src/*.h")

set (Foo_INCLUDE_DIRS "")
foreach (_headerFile ${Foo_HEADERS})
    get_filename_component(_dir ${_headerFile} PATH)
    list (APPEND Foo_INCLUDE_DIRS ${_dir})
endforeach()
list(REMOVE_DUPLICATES Foo_INCLUDE_DIRS)

add_executable (matrixmul ${Foo_SOURCES})

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
