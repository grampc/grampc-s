# This file is part of GRAMPC-S - (https://github.com/grampc/grampc-s)
#
# GRAMPC-S -- A software framework for stochastic model predictive control (SMPC)
#
# Copyright 2024 Daniel Landgraf, Andreas Voelz, Knut Graichen.
# All rights reserved.
#
# GRAMPC-S is distributed under the BSD-3-Clause license, see LICENSE.txt

add_executable(example_double_integrator_GP_prediction
    double_integrator_GP_problem_description.cpp
    double_integrator_GP_prediction_main.cpp
)

add_executable(example_double_integrator_GP_simulation
    double_integrator_GP_problem_description.cpp
    double_integrator_GP_simulation_main.cpp
    true_system_dynamics.cpp
)

target_link_libraries(example_double_integrator_GP_prediction
    grampc-s
)

target_link_libraries(example_double_integrator_GP_simulation
    grampc-s
)
