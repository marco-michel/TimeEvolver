cmake_minimum_required(VERSION 3.13 FATAL_ERROR)

set (CMAKE_CXX_STANDARD 11)

set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_VERBOSE_MAKEFILE ON)


enable_language(CXX)
enable_language(C)

SET(GCC_COVERAGE_COMPILE_FLAGS "-DMKL_ILP64 -m64 -std=c++11 -Wall")

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

#Optional options
#SET(MKLROOT $ENV{MKLROOT})
#SET(CMAKE_EXE_LINKER_FLAGS, "-static")


find_package(Boost COMPONENTS program_options REQUIRED)


#SET HDF5_DIR for MAC usage. Replace the ellipsis with your corresponding prefix directory
#SET(HDF5_DIR .../HDF_Group/HDF5/1.12.0/share/cmake/hdf5)

find_package(HDF5 COMPONENTS CXX REQUIRED)


#USE following 2 lines WITH find_package(MKL)
#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
#include(check_param)


#REMOVE following 2 lines WITH find_package(MKL)
set(BLA_VENDOR Intel10_64ilp)
find_package(BLAS REQUIRED)


#USE following 2 lines WITH find_package(MKL)
#find_package(MKL MODULE REQUIRED)
#set(BLAS_LIBRARIES ${MKL_LINK_PREFIX})


if(DEFINED ENV{MKLROOT})
message("MKLROOT defined")
ELSE()
include_directories(/usr/include/mkl)
ENDIF()

include_directories(${HDF5_INCLUDE_DIR})


set(_hdf5_libs hdf5 hdf5_cpp)
project (TimeEvolver)

add_subdirectory(core TimeEvolver)
add_subdirectory(helper Helper)
add_subdirectory(example Example)




