add_library(TimeEvolver SHARED krylovTimeEvolver.cpp)

target_link_libraries(TimeEvolver PUBLIC
    ${BLAS_LIBRARIES} dl pthread)

target_link_libraries(TimeEvolver PUBLIC  ${CUDA_LIBRARIES})

install(TARGETS TimeEvolver DESTINATION TimeEvolver/lib)
install(FILES matrixDataTypes.h krylovTimeEvolver.h DESTINATION TimeEvolver/include)

if(HDF5_FOUND)
    target_link_libraries(TimeEvolver PUBLIC
    ${HDF5_LIBRARIES}) 
endif()
