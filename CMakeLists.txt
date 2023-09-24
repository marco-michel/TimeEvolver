cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

set (CMAKE_CXX_STANDARD 17)
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



set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/TimeEvolver/lib")
message("Installation directory would be: ${CMAKE_INSTALL_PREFIX}/TimeEvolver/lib")
set (Boost_USE_STATIC_LIBS ON)

include(CTest)

if(NOSYSTEMBOOST)
  set (Boost_NO_SYSTEM_PATHS ON)
endif()

find_package(Boost 1.75 REQUIRED COMPONENTS program_options log log_setup )

find_package(HDF5 COMPONENTS CXX)
if(HDF5_FOUND)
  add_compile_definitions(USE_HDF)
  include_directories(${HDF5_INCLUDE_DIR})
  set(_hdf5_libs hdf5 hdf5_cpp)
endif()


If(APPLE)
  add_compile_definitions(ON_APPLE)
  add_compile_definitions(ACCELERATE_NEW_LAPACK)
  add_compile_definitions(ACCELERATE_LAPACK_ILP64)
  set(BLA_VENDOR Apple)
  find_package(BLAS REQUIRED)
  include_directories(/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks/Accelerate.framework/Headers)
else()
  set(BLA_VENDOR Intel10_64ilp)
  find_package(BLAS)
  if(BLAS_FOUND)
    add_compile_definitions(USE_MKL)
  endif()
  if(NOT BLAS_FOUND)
    unset(BLA_VENDOR)
    find_package(BLAS REQUIRED)
    add_compile_definitions(USE_OPENBLAS)
    set(USING_OPENBLAS TRUE)
    if(USING_ARMADILLO)
      find_package(Armadillo REQUIRED)
      include_directories(${ARMADILLO_INCLUDE_DIRS})
      add_compile_definitions(USE_ARMADILLO)
    endif()
  endif()

endif()


if(DEFINED ENV{MKLROOT})
message("MKLROOT defined")
ELSE()
include_directories(/usr/include/mkl)
ENDIF()

include_directories(${Boost_INCLUDE_DIRS} ${boost_SOURCE_DIR})


add_subdirectory(core TimeEvolver)
add_subdirectory(helper Helper)
add_subdirectory(example Example)

add_subdirectory(test Test)

add_test(NAME BlackHoleTest COMMAND Test/testBlackHole)
add_test(NAME BlackHoleTestLB COMMAND Test/testBlackHoleLB)
add_test(NAME SimpleExampleTest COMMAND Test/testSimpleExample)
