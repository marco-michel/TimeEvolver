


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


#Unlike the tests above this one drives the library directly instead of running
#an example as a subprocess, so it has to link against it.
add_executable(TestErrorHandling testErrorHandling.cpp)
set_target_properties(TestErrorHandling PROPERTIES OUTPUT_NAME "testErrorHandling" INSTALL_RPATH_USE_LINK_PATH TRUE)
target_include_directories(TestErrorHandling PRIVATE ../core)
target_link_libraries(TestErrorHandling PRIVATE TimeEvolver MatrixDataTypes)



