# This file is part of GRAMPC-S - (https://github.com/grampc/grampc-s)
#
# GRAMPC-S -- A software framework for stochastic model predictive control (SMPC)
#
# Copyright 2024 Daniel Landgraf, Andreas Voelz, Knut Graichen.
# All rights reserved.
#
# GRAMPC-S is distributed under the BSD-3-Clause license, see LICENSE.txt

add_executable(example_NLChain_simulation
    NLChain_2_problem_description.cpp
    NLChain_3_problem_description.cpp
    NLChain_4_problem_description.cpp
    NLChain_6_problem_description.cpp
    NLChain_8_problem_description.cpp
    NLChain_10_problem_description.cpp
    NLChain_12_problem_description.cpp
    NLChain_14_problem_description.cpp
    NLChain_simulation_main.cpp
)

target_link_libraries(example_NLChain_simulation
    grampc-s
)
