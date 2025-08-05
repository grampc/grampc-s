# This file is part of GRAMPC-S - (https://github.com/grampc/grampc-s)
#
# GRAMPC-S -- A software framework for stochastic model predictive control (SMPC)
#
# Copyright 2024 Daniel Landgraf, Andreas Voelz, Knut Graichen.
# All rights reserved.
#
# GRAMPC-S is distributed under the BSD-3-Clause license, see LICENSE.txt

add_executable(example_PMSM_prediction
    PMSM_problem_description.cpp
    PMSM_prediction_main.cpp
)

add_executable(example_PMSM_simulation
    PMSM_problem_description.cpp
    PMSM_simulation_main.cpp
)

target_link_libraries(example_PMSM_prediction
    grampc-s
)

target_link_libraries(example_PMSM_simulation
    grampc-s
)
