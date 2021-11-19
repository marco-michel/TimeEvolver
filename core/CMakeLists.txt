add_library(TimeEvolver SHARED krylovTimeEvolver.cpp)

target_link_libraries(TimeEvolver PUBLIC
    ${BLAS_LIBRARIES}
    ${HDF5_LIBRARIES}
    dl pthread) 