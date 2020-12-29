cmake_minimum_required(VERSION 3.13.4)
project(nasieji_2)

set(CMAKE_CXX_STANDARD 17)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/demandPoints.dat
        ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})
include_directories(/usr/include/mpi/)

add_executable(nasieji_2 main.cpp)
target_link_libraries(nasieji_2 ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
    set_target_properties(nasieji_2 PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(nasieji_2 PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
