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

if(TIMEEVOLVER_ENABLE_PYTHON_GTEST)
    if(NOT TIMEEVOLVER_PYTHON)
        message(FATAL_ERROR "TIMEEVOLVER_ENABLE_PYTHON_GTEST requires TIMEEVOLVER_PYTHON=ON.")
    endif()

    add_executable(TestPythonBindingsGTest testPythonBindingsGTest.cpp)
    set_target_properties(TestPythonBindingsGTest PROPERTIES OUTPUT_NAME "testPythonBindingsGTest" INSTALL_RPATH_USE_LINK_PATH TRUE)
    target_include_directories(TestPythonBindingsGTest PRIVATE ../core ../helper)
    target_link_libraries(TestPythonBindingsGTest PRIVATE GTest::gtest GTest::gtest_main TimeEvolver MatrixDataTypes)
    if(HDF5_FOUND)
        target_link_libraries(TestPythonBindingsGTest PRIVATE ${HDF5_LIBRARIES})
    endif()
    target_compile_definitions(TestPythonBindingsGTest PRIVATE TE_PYTHON_MODULE_DIR="$<TARGET_FILE_DIR:timeevolver>")
    add_dependencies(TestPythonBindingsGTest timeevolver)
    gtest_discover_tests(TestPythonBindingsGTest TEST_PREFIX "PythonBindingsGTest.")
endif()
