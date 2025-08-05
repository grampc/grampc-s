# This file is part of GRAMPC-S - (https://github.com/grampc/grampc-s)
#
# GRAMPC-S -- A software framework for stochastic model predictive control (SMPC)
#
# Copyright 2024 Daniel Landgraf, Andreas Voelz, Knut Graichen.
# All rights reserved.
#
# GRAMPC-S is distributed under the BSD-3-Clause license, see LICENSE.txt

target_sources(grampc-s PRIVATE 
	Hermite_polynomial_generator.cpp
	Legendre_polynomial_generator.cpp
	multivariate_polynomial.cpp
	polynomial.cpp)