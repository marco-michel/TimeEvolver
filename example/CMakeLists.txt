


add_executable(Example main.cpp exampleHamiltonian.cpp)
set_target_properties(Example PROPERTIES OUTPUT_NAME "main" INSTALL_RPATH_USE_LINK_PATH TRUE)
target_include_directories(Example PUBLIC ../core ../helper)

target_link_libraries(Example PUBLIC
    ${BLAS_LIBRARIES}
    ${Boost_LIBRARIES}
    dl pthread
    Helper TimeEvolver) 


add_executable(SimpleExample simpleExample.cpp)
set_target_properties(SimpleExample PROPERTIES OUTPUT_NAME "SimpleExample")
target_include_directories(SimpleExample PUBLIC ../core ../helper)

target_link_libraries(SimpleExample PUBLIC
    Helper TimeEvolver) 

install(TARGETS Example RUNTIME DESTINATION TimeEvolver/bin/Example)
install(FILES exampleHamiltonian.h DESTINATION TimeEvolver/include)

if(HDF5_FOUND)
    target_link_libraries(Example PUBLIC
    ${HDF5_LIBRARIES}) 
endif()

