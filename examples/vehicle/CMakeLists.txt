add_executable(example_vehicle
    vehicle_problem_description.cpp
    vehicle_main.cpp
)

target_link_libraries(example_vehicle
    grampc-s
)

