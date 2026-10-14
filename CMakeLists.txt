cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)
enable_language(CXX)
# C is not used for any source file, but find_package(BLAS) performs its
# probe (e.g. "Looking for sgemm_") as a C compile check.
enable_language(C)

option (DOWNLOAD_BOOST "Whether or not Boost should be downloaded" OFF)

set(NOSYSTEMBOOST OFF)


if(DOWNLOAD_BOOST)
  project(EXTERNAL NONE)
  include (cmake/CMakeLists.txt)
  return()
else()
  project(TimeEvolver)
endif()



option (ENABLE_WARNINGS "Compile with the common warning set" ON)
if(ENABLE_WARNINGS AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  # -Wno-unused-parameter: the observable interface deliberately ignores
  # arguments in some overrides.
  add_compile_options(-Wall -Wextra -Wno-unused-parameter)
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
  message(STATUS "HDF5 found - building with HDF5 support (USE_HDF)")
  add_compile_definitions(USE_HDF)
  include_directories(${HDF5_INCLUDE_DIRS})
else()
  message(STATUS "HDF5 not found - HDF5 output and the HDF5 test paths are disabled")
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

# BUILD_TESTING is defined by include(CTest) above and defaults to ON.
# The individual tests are registered by tests/CMakeLists.txt.
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
