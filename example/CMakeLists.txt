


add_executable(Example main.cpp exampleHamiltonian.cpp)
set_target_properties(Example PROPERTIES OUTPUT_NAME "main")
target_include_directories(Example PUBLIC ../core ../helper)

target_link_libraries(Example PUBLIC
    ${BLAS_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${Boost_LIBRARIES}
    dl pthread
    Helper TimeEvolver) 

add_executable(SimpleExample simpleExample.cpp)
set_target_properties(SimpleExample PROPERTIES OUTPUT_NAME "SimpleExample")
target_include_directories(SimpleExample PUBLIC ../core ../helper)

target_link_libraries(SimpleExample PUBLIC
    Helper TimeEvolver) 