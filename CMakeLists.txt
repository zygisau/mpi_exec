cmake_minimum_required(VERSION 3.13.4)
project(nasieji_2)

set(CMAKE_CXX_STANDARD 17)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/demandPoints.dat
        ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
#include_directories(/usr/include/mpi/)

SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)

set(SOURCE_FILES main.cpp)
add_executable(nasieji_2 ${SOURCE_FILES})
