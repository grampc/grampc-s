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


#ifndef HERMITE_POLYNOMIAL_GENERATOR_HPP
#define HERMITE_POLYNOMIAL_GENERATOR_HPP

#include <algorithm>
#include "orthogonal_polynomial_generator.hpp"
#include "util/grampc_s_util.hpp"

namespace grampc
{
    // Generator for univariate Hermite polynomials with measure 1/sqrt(2*pi) * exp(-x^2 / 2) 
    class HermitePolynomialGenerator : public OrthogonalPolynomialGenerator
    {
    public:
        // Constructor specifying the maximum polynomial order
        HermitePolynomialGenerator(typeInt maxOrder);

        // Get a polynomial
        virtual PolynomialConstPtr getPolynomial(typeInt order) const override;

        // Get the squared norm of the polynomial 
        virtual typeRNum getSquaredNorm(typeInt order) const override;

        // Get maximum order of polynomial that can be created
        virtual typeInt getMaximumOrder() const override;

    private:
        std::vector<PolynomialPtr> hermitePolynomials;
        Vector hermiteSquaredNorm;
    };   
}


#endif //HERMITE_POLYNOMIAL_GENERATOR_HPP