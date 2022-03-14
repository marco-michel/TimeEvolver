cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

set (CMAKE_CXX_STANDARD 11)
enable_language(CXX)
enable_language(C)

set(CMAKE_VERBOSE_MAKEFILE ON)

option (DOWNLOAD_BOOST "Whether or not Boost should be downloaded" OFF)

set(NOSYSTEMBOOST OFF)


if(DOWNLOAD_BOOST)
  project(EXTERNAL NONE)
  include (cmake/CMakeLists.txt)
  return()
else()
  project(TimeEvolver)
endif()


set (Boost_USE_STATIC_LIBS ON)

if(NOSYSTEMBOOST)
  set (Boost_NO_SYSTEM_PATHS ON)
endif()

find_package(Boost 1.75 REQUIRED COMPONENTS program_options)

find_package(HDF5 COMPONENTS CXX REQUIRED)

set(BLA_VENDOR Intel10_64ilp)
find_package(BLAS REQUIRED)

if(DEFINED ENV{MKLROOT})
message("MKLROOT defined")
ELSE()
include_directories(/usr/include/mkl)
ENDIF()

include_directories(${HDF5_INCLUDE_DIR} ${Boost_INCLUDE_DIRS} ${boost_SOURCE_DIR})


set(_hdf5_libs hdf5 hdf5_cpp)


add_subdirectory(core TimeEvolver)
add_subdirectory(helper Helper)
add_subdirectory(example Example)
