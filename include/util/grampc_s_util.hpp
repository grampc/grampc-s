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

#ifndef GRAMPC_S_UTIL_HPP
#define GRAMPC_S_UTIL_HPP

#include <iostream>
#include <fstream>
#include "grampc_s_constants.hpp"
#include "point_transformation/quadrature_rules/Hermite_quadrature.hpp"
#include "point_transformation/quadrature_rules/Legendre_quadrature.hpp"
#include "polynomial/Hermite_polynomial_generator.hpp"
#include "polynomial/Legendre_polynomial_generator.hpp"
#include "grampc_interface/grampc_interface.hpp"

namespace grampc
{
    // Factorial
    long long factorial(long long n);

    // Return the quadrature rule that corresponds to the distribution
    QuadratureRuleConstPtr correspondingQuadratureRule(PolynomialFamily polyFam, typeInt quadratureOrder);

    // Return the quadrature rule that corresponds to the distribution
    PolynomialGeneratorPtr correspondingPolynomialGenerator(PolynomialFamily polyFam, typeInt maxOrder);

    // Derivative of the Cholesky decomposition
    void deriveCholesky(VectorRef out, MatrixConstRef chol, typeInt k, typeInt l);

    // Writes the solution of the solver in text files
    void writeTrajectoriesToFile(GrampcPtr solver, typeInt numberOfStates);
  }

#endif // GRAMPC_S_UTIL_HPP
