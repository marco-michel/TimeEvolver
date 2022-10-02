


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



