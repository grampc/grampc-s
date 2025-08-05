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


#ifndef POLYNOMIAL_GENERATOR_HPP
#define POLYNOMIAL_GENERATOR_HPP

#include "polynomial.hpp"

namespace grampc
{
    // Generator for a family of orthogonal univariate polynomials
    class OrthogonalPolynomialGenerator
    {
    public:
        virtual ~OrthogonalPolynomialGenerator() {};

        // Get a polynomial
        virtual PolynomialConstPtr getPolynomial(typeInt order) const = 0;

        // Get the squared norm of the polynomial 
        virtual typeRNum getSquaredNorm(typeInt order) const = 0;

        // Get maximum order of polynomial that can be created
        virtual typeInt getMaximumOrder() const = 0;
    };   

    // Alias
    typedef std::shared_ptr<OrthogonalPolynomialGenerator> PolynomialGeneratorPtr;
}

#endif // POLYNOMIAL_GENERATOR_HPP