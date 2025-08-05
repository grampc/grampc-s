# This file is part of GRAMPC-S - (https://github.com/grampc/grampc-s)
#
# GRAMPC-S -- A software framework for stochastic model predictive control (SMPC)
#
# Copyright 2024 Daniel Landgraf, Andreas Voelz, Knut Graichen.
# All rights reserved.
#
# GRAMPC-S is distributed under the BSD-3-Clause license, see LICENSE.txt

add_executable(example_inverted_pendulum_prediction
    inverted_pendulum_problem_description.cpp
    inverted_pendulum_prediction_main.cpp
)

add_executable(example_inverted_pendulum_simulation
    inverted_pendulum_problem_description.cpp
    inverted_pendulum_simulation_main.cpp
)

target_link_libraries(example_inverted_pendulum_prediction
   grampc-s
)

target_link_libraries(example_inverted_pendulum_simulation
    grampc-s
)