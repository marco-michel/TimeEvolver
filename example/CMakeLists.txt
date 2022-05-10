


add_executable(Example main.cpp exampleHamiltonian.cpp)
set_target_properties(Example PROPERTIES OUTPUT_NAME "main")
target_include_directories(Example PUBLIC ../core ../helper)

target_link_libraries(Example PUBLIC
    ${BLAS_LIBRARIES}
    ${Boost_LIBRARIES}
    dl pthread
    Helper TimeEvolver) 

if(HDF5_FOUND)
    target_link_libraries(Example PUBLIC
    ${HDF5_LIBRARIES}) 
endif()