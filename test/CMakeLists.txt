


add_executable(Test test.cpp)
set_target_properties(Test PROPERTIES OUTPUT_NAME "test" INSTALL_RPATH_USE_LINK_PATH TRUE)

if(HDF5_FOUND)
    target_link_libraries(Test PUBLIC
    ${HDF5_LIBRARIES}) 
endif()



