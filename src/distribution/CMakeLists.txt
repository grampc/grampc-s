# This file is part of GRAMPC-S - (https://github.com/grampc/grampc-s)
#
# GRAMPC-S -- A software framework for stochastic model predictive control (SMPC)
#
# Copyright 2024 Daniel Landgraf, Andreas Voelz, Knut Graichen.
# All rights reserved.
#
# GRAMPC-S is distributed under the BSD-3-Clause license, see LICENSE.txt

target_sources(grampc-s PRIVATE 
	distribution.cpp
	univariate_uniform_distribution.cpp
	univariate_piecewise_constant_distribution.cpp
	multivariate_uncorrelated_distribution.cpp
	Gaussian_distribution.cpp
	univariate_gamma_distribution.cpp
	univariate_lognormal_distribution.cpp
	univariate_student-t_distribution.cpp
	univariate_weibull_distribution.cpp
	univariate_exponential_distribution.cpp
	univariate_extreme_value_distribution.cpp
	univariate_chi_squared_distribution.cpp
	univariate_f_distribution.cpp
	univariate_beta_distribution.cpp)