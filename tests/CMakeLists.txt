include(FetchContent)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.zip
)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)
include(GoogleTest)


add_executable(unit_tests test_smatrix_equals.cpp test_smatrix_creation.cpp ../example/exampleHamiltonian.cpp)
target_include_directories(unit_tests PRIVATE ../core ../helper ../example)
target_link_libraries(unit_tests PRIVATE GTest::gtest_main MatrixDataTypes Helper ${HDF5_LIBRARIES})



add_executable(TestBlackHole testBlackHole.cpp)
set_target_properties(TestBlackHole PROPERTIES OUTPUT_NAME "testBlackHole" INSTALL_RPATH_USE_LINK_PATH TRUE)

if(HDF5_FOUND)
    target_link_libraries(TestBlackHole PUBLIC
    ${HDF5_LIBRARIES}) 
endif()


add_executable(TestBlackHoleLB testBlackHoleLB.cpp)
set_target_properties(TestBlackHoleLB PROPERTIES OUTPUT_NAME "testBlackHoleLB" INSTALL_RPATH_USE_LINK_PATH TRUE)

if(HDF5_FOUND)
    target_link_libraries(TestBlackHoleLB PUBLIC
    ${HDF5_LIBRARIES}) 
endif()


add_executable(TestSimpleExample testSimpleExample.cpp)
set_target_properties(TestSimpleExample PROPERTIES OUTPUT_NAME "testSimpleExample" INSTALL_RPATH_USE_LINK_PATH TRUE)



add_test(NAME BlackHoleTest COMMAND testBlackHole)
add_test(NAME BlackHoleTestLB COMMAND testBlackHoleLB)
add_test(NAME SimpleExampleTest COMMAND testSimpleExample)
gtest_discover_tests(unit_tests)