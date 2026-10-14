pybind11_add_module(timeevolver bindings.cpp)

target_link_libraries(timeevolver PRIVATE TimeEvolver)
target_include_directories(timeevolver PRIVATE
    ${CMAKE_SOURCE_DIR}/core
    ${CMAKE_SOURCE_DIR}/helper
)
