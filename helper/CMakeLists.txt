
add_library(Helper STATIC Basis.cpp hamiltonian.cpp)
add_library(outputHelper STATIC krylovHelper.cpp)

target_include_directories(Helper PUBLIC ../core)
target_include_directories(outputHelper PUBLIC ../core)

set_target_properties(Helper PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

target_link_libraries(Helper PUBLIC
    ${Boost_LIBRARIES}
    dl pthread MatrixDataTypes)  


target_link_libraries(outputHelper PUBLIC
    dl pthread MatrixDataTypes)  

if(HDF5_FOUND)
    target_link_libraries(outputHelper PUBLIC
    ${HDF5_LIBRARIES}) 
endif()    

install(TARGETS Helper LIBRARY DESTINATION TimeEvolver/lib)    
install(TARGETS outputHelper LIBRARY DESTINATION TimeEvolver/lib)    
install(FILES Basis.h hamiltonian.h krylovHelper.h parameter.h DESTINATION TimeEvolver/include)

