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
#include "point_transformation/quadrature_rules/Legendre_quadrature.hpp"

namespace grampc
{
    LegendreQuadrature::LegendreQuadrature(typeInt quadratureOrder)
      : pointsNormalized_(quadratureOrder),
        roots_(quadratureOrder),
        weights_(quadratureOrder),
        quadratureOrder_(quadratureOrder)
    {
        LegendrePolynomialGenerator polyGen = LegendrePolynomialGenerator(quadratureOrder);
        Eigen::PolynomialSolver<typeRNum, Eigen::Dynamic> polySolver;

        // Compute Legendre polynomial 
        PolynomialConstPtr poly = polyGen.getPolynomial(quadratureOrder);

        // Compute roots of the Legendre polynomial
        polySolver.compute(poly->getCoefficients());
        roots_ = polySolver.roots().real();

        // Normalized points for a univariate distriubtion with mean = 0 and variance = 1
        pointsNormalized_ = std::sqrt(3.0) * roots_;
        
        // Compute weights for the univariate quadrature points
        // Measure for Legendre: 1 / 2
        for (int j = 0; j < quadratureOrder; ++j)
        {
            weights_[j] = 1.0 / ((1.0 - std::pow(roots_[j], 2)) * std::pow(polyGen.getPolynomial(quadratureOrder_)->gradient(roots_[j]), 2));
        }
    }

    const Vector& LegendreQuadrature::polynomialRoots() const
    {
        return roots_;
    }

    const Vector& LegendreQuadrature::pointsNormalized() const
    {
        return pointsNormalized_;
    }

    const Vector& LegendreQuadrature::weights() const
    {
        return weights_;
    }

    typeInt LegendreQuadrature::numberOfPoints() const
    {
        return quadratureOrder_;
    }
}