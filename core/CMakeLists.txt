
add_library(MatrixDataTypes STATIC matrixDataTypes.cpp)
target_link_libraries(MatrixDataTypes PUBLIC ${BLAS_LIBRARIES} ${Boost_LIBRARIES}
     dl pthread)


add_library(TimeEvolver STATIC krylovTimeEvolver.cpp krylovObservables.cpp krylovLogger.cpp)
target_include_directories(TimeEvolver PUBLIC ../helper)

target_link_libraries(TimeEvolver PUBLIC ${BLAS_LIBRARIES} ${Boost_LIBRARIES}
     dl pthread MatrixDataTypes)

install(TARGETS TimeEvolver DESTINATION TimeEvolver/lib)
install(FILES matrixDataTypes.h krylovTimeEvolver.h krylovObservables.h mathHeader.h DESTINATION TimeEvolver/include)

if(HDF5_FOUND)
    target_link_libraries(TimeEvolver PUBLIC
    ${HDF5_LIBRARIES}) 
endif()

if(USING_OPENBLAS)
    target_link_libraries(TimeEvolver PUBLIC lapacke) 
endif()

if(USING_ARMADILLO)
target_link_libraries(TimeEvolver PUBLIC  ${ARMADILLO_LIBRARIES})
endif()