add_executable(TestBlackHole testBlackHole.cpp)
set_target_properties(TestBlackHole PROPERTIES OUTPUT_NAME "testBlackHole" INSTALL_RPATH_USE_LINK_PATH TRUE)
target_link_libraries(TestBlackHole PRIVATE GTest::gtest_main)
if(HDF5_FOUND)
    target_link_libraries(TestBlackHole PRIVATE ${HDF5_LIBRARIES})
endif()
gtest_discover_tests(
    TestBlackHole
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_executable(TestBlackHoleLB testBlackHoleLB.cpp)
set_target_properties(TestBlackHoleLB PROPERTIES OUTPUT_NAME "testBlackHoleLB" INSTALL_RPATH_USE_LINK_PATH TRUE)
target_link_libraries(TestBlackHoleLB PRIVATE GTest::gtest_main)
if(HDF5_FOUND)
    target_link_libraries(TestBlackHoleLB PRIVATE ${HDF5_LIBRARIES})
endif()
gtest_discover_tests(
    TestBlackHoleLB
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_executable(TestSimpleExample testSimpleExample.cpp)
set_target_properties(TestSimpleExample PROPERTIES OUTPUT_NAME "testSimpleExample" INSTALL_RPATH_USE_LINK_PATH TRUE)
target_link_libraries(TestSimpleExample PRIVATE GTest::gtest_main)
gtest_discover_tests(
    TestSimpleExample
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

add_executable(TestBackendParity testBackendParity.cpp ../example/exampleHamiltonian.cpp)
set_target_properties(TestBackendParity PROPERTIES OUTPUT_NAME "testBackendParity" INSTALL_RPATH_USE_LINK_PATH TRUE)
target_include_directories(TestBackendParity PUBLIC ../core ../helper ../example)
target_link_libraries(TestBackendParity PRIVATE GTest::gtest_main Helper TimeEvolver)
gtest_discover_tests(
    TestBackendParity
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    PROPERTIES SKIP_REGULAR_EXPRESSION ".*SKIPPED.*"
)
