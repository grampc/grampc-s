# This file is part of GRAMPC-S - (https://github.com/grampc/grampc-s)
#
# GRAMPC-S -- A software framework for stochastic model predictive control (SMPC)
#
# Copyright 2024 Daniel Landgraf, Andreas Voelz, Knut Graichen.
# All rights reserved.
#
# GRAMPC-S is distributed under the BSD-3-Clause license, see LICENSE.txt

target_sources(grampc-s PRIVATE 
	monte_carlo_problem_description.cpp
	sigma_point_problem_description.cpp
	problem_description.cpp
	taylor_problem_description.cpp
	resampling_problem_description.cpp
	resampling_GP_problem_description.cpp)