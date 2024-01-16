add_executable(example_double_integrator
    double_integrator_problem_description.cpp
    double_integrator_main.cpp
)

target_link_libraries(example_double_integrator
    grampc-s
)
