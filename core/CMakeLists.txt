add_library(TimeEvolver SHARED krylovTimeEvolver.cpp krylovObservables.cpp krylovLogger.cpp matrixDataTypes.cpp)
target_include_directories(TimeEvolver PUBLIC ../helper)

target_link_libraries(TimeEvolver PUBLIC
    ${BLAS_LIBRARIES} dl pthread)

install(TARGETS TimeEvolver DESTINATION TimeEvolver/lib)
install(FILES matrixDataTypes.h krylovTimeEvolver.h krylovObservables.h mathHeader.h DESTINATION TimeEvolver/include)

if(HDF5_FOUND)
    target_link_libraries(TimeEvolver PUBLIC
    ${HDF5_LIBRARIES}) 
endif()
