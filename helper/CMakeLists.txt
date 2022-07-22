
add_library(Helper SHARED Basis.cpp hamiltonian.cpp)
target_include_directories(Helper PUBLIC ../core)

set_target_properties(Helper PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

target_link_libraries(Helper PUBLIC
    ${BLAS_LIBRARIES}
    ${Boost_LIBRARIES}
    dl pthread) 

if(OpenMP_CXX_FOUND)
    target_link_libraries(Helper PUBLIC OpenMP::OpenMP_CXX)
endif()

if(HDF5_FOUND)
    target_link_libraries(Helper PUBLIC
    ${HDF5_LIBRARIES}) 
endif()    

install(TARGETS Helper LIBRARY DESTINATION TimeEvolver/lib)    
install(FILES Basis.h hamiltonian.h DESTINATION TimeEvolver/include)

