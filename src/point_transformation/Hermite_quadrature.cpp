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


#include <unsupported/Eigen/Polynomials>
#include "point_transformation/quadrature_rules//Hermite_quadrature.hpp"

namespace grampc
{
    HermiteQuadrature::HermiteQuadrature(typeInt quadratureOrder)
      : roots_(quadratureOrder),
        weights_(quadratureOrder),
        quadratureOrder_(quadratureOrder)
    {
        HermitePolynomialGenerator polyGen = HermitePolynomialGenerator(quadratureOrder);
        Eigen::PolynomialSolver<typeRNum, Eigen::Dynamic> polySolver;

        // Compute Hermite polynomial 
        PolynomialConstPtr poly = polyGen.getPolynomial(quadratureOrder);

        // Compute roots of the Hermite polynomial
        polySolver.compute(poly->getCoefficients());
        roots_ = polySolver.roots().real();
        
        // Compute weights for the univariate quadrature points
        // Measure for Hermite: 1/sqrt(2*pi) * exp(-x^2 / 2)
        for (int j = 0; j < quadratureOrder; ++j)
        {
            weights_[j] = ((typeRNum)factorial(quadratureOrder_)) / (std::pow(quadratureOrder_ * polyGen.getPolynomial(quadratureOrder - 1)->evaluate(roots_[j]), 2));
        }
    }

    const Vector& HermiteQuadrature::polynomialRoots() const
    {
        return roots_;
    }

    const Vector& HermiteQuadrature::pointsNormalized() const
    {
        return roots_;
    }

    const Vector& HermiteQuadrature::weights() const
    {
        return weights_;
    }

    typeInt HermiteQuadrature::numberOfPoints() const
    {
        return quadratureOrder_;
    }
}