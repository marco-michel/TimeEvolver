
add_executable(Benchmark benchmark.cpp ../example/exampleHamiltonian.cpp)
set_target_properties(Benchmark PROPERTIES OUTPUT_NAME "benchmark" INSTALL_RPATH_USE_LINK_PATH TRUE)
target_include_directories(Benchmark PUBLIC ../core ../helper ../example)

target_link_libraries(Benchmark PUBLIC
    ${BLAS_LIBRARIES}
    ${Boost_LIBRARIES}
    dl pthread
    Helper TimeEvolver)

if(HDF5_FOUND)
    target_link_libraries(Benchmark PUBLIC ${HDF5_LIBRARIES})
endif()

install(TARGETS Benchmark RUNTIME DESTINATION TimeEvolver/bin/Benchmark)
