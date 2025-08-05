# This file is part of GRAMPC-S - (https://github.com/grampc/grampc-s)
#
# GRAMPC-S -- A software framework for stochastic model predictive control (SMPC)
#
# Copyright 2024 Daniel Landgraf, Andreas Voelz, Knut Graichen.
# All rights reserved.
#
# GRAMPC-S is distributed under the BSD-3-Clause license, see LICENSE.txt

target_sources(grampc-s PRIVATE 
	composed_quadrature.cpp
	Hermite_quadrature.cpp
	Legendre_quadrature.cpp
	monte_carlo.cpp
	PCE_transformation.cpp
	stirling_interpolation_first_order.cpp
	stirling_interpolation_second_order.cpp
	unscented_transformation.cpp)