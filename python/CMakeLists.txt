find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 REQUIRED)

pybind11_add_module(grampc_s MODULE)
target_link_libraries(grampc_s PRIVATE Eigen3::Eigen grampc-s)

include_directories(include)

add_subdirectory(src)

install(TARGETS grampc_s DESTINATION .)