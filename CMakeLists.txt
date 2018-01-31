cmake_minimum_required(VERSION 2.8)

project(SPIN_BOSON Fortran)

#set(CMAKE_Fortran_COMPILER      "mpif90")
#set(CMAKE_Fortran_FLAGS   "-O2")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR})

set(TARGET_NAME   "shin-metiu-model")

set(GLOBAL_LIB "global_lib")
set(TDSE_LIB   "tdse_lib")
set(BO_LIB   "bo_lib")
set(IO_LIB   "io_lib")

set(SRC_DIR "./src")
add_subdirectory(${SRC_DIR})
