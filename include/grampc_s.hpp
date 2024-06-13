/* This file is part of GRAMPC-S - (https://github.com/grampc/grampc-s)
 *
 * GRAMPC-S -- A software framework for stochastic model predictive control (SMPC)
 *
 * Copyright 2024 Daniel Landgraf, Andreas Voelz, Knut Graichen.
 * All rights reserved.
 *
 * GRAMPC-S is distributed under the BSD-3-Clause license, see LICENSE.txt
 *
 */


#ifndef GRAMPC_S_HPP
#define GRAMPC_S_HPP

// headers from the c++ standard library
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <Eigen/Dense>

// include grampc library
#include "grampc_interface/grampc_interface.hpp"

// include problem descriptions
#include "problem_description/sigma_point_problem_description.hpp"
#include "problem_description/taylor_problem_description.hpp"
#include "problem_description/monte_carlo_problem_description.hpp"
#include "problem_description/resampling_problem_description.hpp"
#include "problem_description/resampling_GP_problem_description.hpp"

// include distributions
#include "distribution/Gaussian_distribution.hpp"
#include "distribution/univariate_chi_squared_distribution.hpp"
#include "distribution/univariate_exponential_distribution.hpp"
#include "distribution/multivariate_uncorrelated_distribution.hpp"
#include "distribution/univariate_piecewise_constant_distribution.hpp"
#include "distribution/univariate_uniform_distribution.hpp"
#include "distribution/univariate_gamma_distribution.hpp"
#include "distribution/univariate_lognormal_distribution.hpp"
#include "distribution/univariate_student-t_distribution.hpp"
#include "distribution/univariate_weibull_distribution.hpp"
#include "distribution/univariate_extreme_value_distribution.hpp"
#include "distribution/univariate_f_distribution.hpp"
#include "distribution/univariate_beta_distribution.hpp"

// incluede polynomials
#include "polynomial/multivariate_polynomial.hpp"
#include "polynomial/orthogonal_polynomial_generator.hpp"
#include "polynomial/Hermite_polynomial_generator.hpp"
#include "polynomial/Legendre_polynomial_generator.hpp"

// inlcude point generators
#include "point_transformation/monte_carlo.hpp"
#include "point_transformation/unscented_transformation.hpp"
#include "point_transformation/stirling_interpolation_first_order.hpp"
#include "point_transformation/stirling_interpolation_second_order.hpp"
#include "point_transformation/quadrature_rules/Legendre_quadrature.hpp"
#include "point_transformation/quadrature_rules/Hermite_quadrature.hpp"
#include "point_transformation/composed_quadrature.hpp"
#include "point_transformation/PCE_transformation.hpp"


// include chance constraint approximations
#include "constraint_approx/chebyshev_constraint_approximation.hpp"
#include "constraint_approx/Gaussian_constraint_approximation.hpp"
#include "constraint_approx/symmetric_constraint_approximation.hpp"

// include Gaussian processes and kernels
#include "gaussian_process/stationary_kernel.hpp"
#include "gaussian_process/gaussian_process.hpp"
#include "gaussian_process/squared_exponential_kernel.hpp"
#include "gaussian_process/periodic_kernel.hpp"
#include "gaussian_process/locally_periodic_kernel.hpp"
#include "gaussian_process/Matern32_kernel.hpp"
#include "gaussian_process/Matern52_kernel.hpp"
#include "gaussian_process/kernel_sum.hpp"
#include "gaussian_process/kernel_product.hpp"

// include simulator
#include "simulator/simulator.hpp"

// include utility functions
#include "util/grampc_s_constants.hpp"
#include "util/grampc_s_util.hpp"

#endif // GRAMPC_S_HPP