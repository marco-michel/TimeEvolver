# GoogleTest has to be provided by the system; it is deliberately not downloaded
# during configuration so that a build never depends on network access.
find_package(GTest QUIET)

# Default the unit tests to whatever GoogleTest availability allows, but fail
# loudly if they were asked for explicitly and cannot be built.
option(TIMEEVOLVER_UNIT_TESTS "Build the GoogleTest based unit tests (requires an installed GoogleTest)" ${GTest_FOUND})

if(TIMEEVOLVER_UNIT_TESTS AND NOT GTest_FOUND)
  message(FATAL_ERROR
    "TIMEEVOLVER_UNIT_TESTS is ON but GoogleTest was not found.\n"
    "Install it (Debian/Ubuntu: 'sudo apt install libgtest-dev', "
    "macOS: 'brew install googletest') or reconfigure with "
    "-DTIMEEVOLVER_UNIT_TESTS=OFF.")
endif()

if(TIMEEVOLVER_UNIT_TESTS)
  include(GoogleTest)
  message(STATUS "GoogleTest found - building the unit tests")
else()
  message(STATUS "GoogleTest not found - the unit tests are disabled")
endif()


if(TIMEEVOLVER_UNIT_TESTS)
  add_executable(unit_tests test_smatrix_equals.cpp test_smatrix_creation.cpp ../example/exampleHamiltonian.cpp)
  target_include_directories(unit_tests PRIVATE ../core ../helper ../example)
  target_link_libraries(unit_tests PRIVATE GTest::gtest_main MatrixDataTypes Helper ${HDF5_LIBRARIES})
endif()



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

if(TIMEEVOLVER_UNIT_TESTS)
  gtest_discover_tests(unit_tests)
endif()

if(TIMEEVOLVER_ENABLE_PYTHON_GTEST)
    if(NOT TIMEEVOLVER_PYTHON)
        message(FATAL_ERROR "TIMEEVOLVER_ENABLE_PYTHON_GTEST requires TIMEEVOLVER_PYTHON=ON.")
    endif()
    if(NOT TIMEEVOLVER_UNIT_TESTS)
        message(FATAL_ERROR "TIMEEVOLVER_ENABLE_PYTHON_GTEST requires TIMEEVOLVER_UNIT_TESTS=ON, i.e. an installed GoogleTest.")
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
