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


#include "polynomial/Legendre_polynomial_generator.hpp"

namespace grampc
{
    LegendrePolynomialGenerator::LegendrePolynomialGenerator(typeInt maxOrder)
    : legendrePolynomials(std::vector<grampc::PolynomialPtr>(std::max(1, maxOrder) + 1)),
      legendreSquaredNorm(Vector::Zero(std::max(1, maxOrder) + 1))
    {
        // vector of polynomial coefficients
        Vector coeffs(2);
        
        /***** Compute coefficients for the Legendre polynomials *****/
        // Legendre polynomial of zero order
        legendrePolynomials[0] = Poly(Vector::Constant(1,1));

        // Legendre polynomial of first order
        coeffs[0] = 0.0;
        coeffs[1] = 1.0;
        legendrePolynomials[1] = Poly(coeffs);

        // number of coefficients
        typeInt numCoeffs = 2;

        // Recursive computation of higher order polynomials
        for(int i = 2; i < maxOrder + 1; ++i)
        {
            coeffs.conservativeResize(++numCoeffs);

            // set coefficients of order zero
            coeffs[0] = - (i - 1.0) / i * legendrePolynomials[i - 2]->getCoefficient(0);

            // compute coefficients of higher order
            for(int j = 2 - i % 2 ; j < numCoeffs - 1; j = j + 2)
            {
                coeffs[j] = - (i - 1.0) / i * legendrePolynomials[i - 2]->getCoefficient(j) + (2.0 * i - 1.0) / i * legendrePolynomials[i - 1]->getCoefficient(j-1);
            }
            for(int j = 1 - i % 2 ; j < numCoeffs - 1; j = j + 2)
            {
                coeffs[j] = 0.0;
            }

            // compute last coefficient
            coeffs[numCoeffs - 1] = (2.0 * i - 1.0) / i * legendrePolynomials[i - 1]->getCoefficient(numCoeffs - 2);

            // set coefficients of the i-th order Legendre polynomial
            legendrePolynomials[i] = Poly(coeffs);
        }


        /***** Compute norms of the polynomials *****/
        for(int i = 0; i < maxOrder + 1; ++ i)
        {
            legendreSquaredNorm[i] = 1.0 / (2.0 * i + 1.0);
        }
    }


    PolynomialConstPtr LegendrePolynomialGenerator::getPolynomial(typeInt order) const
    {
        return legendrePolynomials[order];
    }


    typeRNum LegendrePolynomialGenerator::getSquaredNorm(typeInt order) const
    {
        return legendreSquaredNorm[order];
    }

    typeInt LegendrePolynomialGenerator::getMaximumOrder() const
    {
        return legendrePolynomials.size() - 1;
    }
}