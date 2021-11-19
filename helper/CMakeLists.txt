
add_library(Helper SHARED Basis.cpp hamiltonian.cpp)
target_include_directories(Helper PUBLIC ../core)

target_link_libraries(Helper PUBLIC
    ${BLAS_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${Boost_LIBRARIES}
    dl pthread) 