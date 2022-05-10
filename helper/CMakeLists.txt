
add_library(Helper SHARED Basis.cpp hamiltonian.cpp)
target_include_directories(Helper PUBLIC ../core)

target_link_libraries(Helper PUBLIC
    ${BLAS_LIBRARIES}
    ${Boost_LIBRARIES}
    dl pthread) 

if(HDF5_FOUND)
    target_link_libraries(Helper PUBLIC
    ${HDF5_LIBRARIES}) 
endif()