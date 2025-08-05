# This file is part of GRAMPC-S - (https://github.com/grampc/grampc-s)
#
# GRAMPC-S -- A software framework for stochastic model predictive control (SMPC)
#
# Copyright 2024 Daniel Landgraf, Andreas Voelz, Knut Graichen.
# All rights reserved.
#
# GRAMPC-S is distributed under the BSD-3-Clause license, see LICENSE.txt

target_sources(grampc-s PRIVATE 
	squared_exponential_kernel.cpp
	gaussian_process.cpp
	periodic_kernel.cpp
	locally_periodic_kernel.cpp
	Matern32_kernel.cpp
	Matern52_kernel.cpp
	kernel_sum.cpp
	kernel_product.cpp)