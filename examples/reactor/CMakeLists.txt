add_executable(example_reactor
    reactor_problem_description.cpp
    reactor_main.cpp
)

target_link_libraries(example_reactor
    grampc-s
)
