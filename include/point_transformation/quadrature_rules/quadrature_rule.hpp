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

#ifndef QUADRATURE_RULE_HPP
#define QUADRATURE_RULE_HPP

#include "util/grampc_s_constants.hpp"


namespace grampc
{
    // Univariate quadrature rule
    class QuadratureRule
    {
    public:
        virtual ~QuadratureRule() {};

        // Roots of the corresponding orthogonal polynomials
        virtual const Vector& polynomialRoots() const = 0;

        // Get points for a normalized distribution with zero mean and variance = 1
        virtual const Vector& pointsNormalized() const = 0;
        
        // Get weights of the quadrature points 
        virtual const Vector& weights() const = 0;

        // Number of quadrature points
        virtual typeInt numberOfPoints() const = 0;
    };   

    // Alias
    typedef std::shared_ptr<QuadratureRule> QuadratureRulePtr;
    typedef std::shared_ptr<const QuadratureRule> QuadratureRuleConstPtr;
}

#endif // QUADRATURE_RULE_HPP