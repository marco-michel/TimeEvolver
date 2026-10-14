
add_library(MatrixDataTypes STATIC matrixDataTypes.cpp)
target_link_libraries(MatrixDataTypes PUBLIC ${BLAS_LIBRARIES} ${Boost_LIBRARIES}
     dl pthread)

if(TIMEEVOLVER_PYTHON)
    set_target_properties(MatrixDataTypes PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()


add_library(TimeEvolver STATIC krylovTimeEvolver.cpp krylovObservables.cpp krylovLogger.cpp)
target_include_directories(TimeEvolver PUBLIC ../helper)

target_link_libraries(TimeEvolver PUBLIC ${BLAS_LIBRARIES} ${Boost_LIBRARIES}
     dl pthread MatrixDataTypes)

if(TIMEEVOLVER_PYTHON)
    set_target_properties(TimeEvolver PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

install(TARGETS TimeEvolver DESTINATION TimeEvolver/lib)
install(FILES matrixDataTypes.h krylovTimeEvolver.h krylovObservables.h krylovSampleWriter.h mathHeader.h krylovExceptions.h hdf5Support.h version.h DESTINATION TimeEvolver/include)

if(USING_OPENBLAS)
    target_link_libraries(TimeEvolver PUBLIC lapacke) 
endif()

if(USING_ARMADILLO)
target_link_libraries(TimeEvolver PUBLIC  ${ARMADILLO_LIBRARIES})
endif()

if(HDF5_FOUND)
    #MatrixDataTypes calls into HDF5 itself, so it has to link it in its own
    #right rather than rely on TimeEvolver pulling the libraries in.
    target_link_libraries(MatrixDataTypes PUBLIC ${HDF5_LIBRARIES})
    target_link_libraries(TimeEvolver PUBLIC ${HDF5_LIBRARIES})
endif()
