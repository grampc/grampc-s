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


#include "polynomial/Hermite_polynomial_generator.hpp"

namespace grampc
{
    HermitePolynomialGenerator::HermitePolynomialGenerator(typeInt maxOrder)
    : hermitePolynomials(std::vector<grampc::PolynomialPtr>(std::max(1, maxOrder) + 1)),
      hermiteSquaredNorm(Vector::Zero(std::max(1, maxOrder) + 1))
    {
        // vector of polynomial coefficients
        Vector coeffs(2);

        /***** Compute coefficients for the probabilist's Hermite polynomials *****/
        // Measure: 1/sqrt(2*pi) * exp(-x^2 / 2)

        // Hermite polynomial of order zero
        hermitePolynomials[0] = Poly(Vector::Constant(1,1));

        // Hermite polynomial of first order
        coeffs[0] = 0.0;
        coeffs[1] = 1.0;
        hermitePolynomials[1] = Poly(coeffs);

        // number of coefficients
        typeInt numCoeffs = 2;

        // Recursive computation of higher order polynomials
        for(int i = 2; i < maxOrder + 1; ++i)
        {
            coeffs.conservativeResize(++numCoeffs);

            // set coefficients of order zero
            coeffs[0] = -(i-1) * hermitePolynomials[i - 2]->getCoefficient(0);

            // compute coefficients of higher order
            for(int j = 1; j < numCoeffs - 2; ++j)
            {
                coeffs[j] = hermitePolynomials[i-1]->getCoefficient(j-1) - (i-1) * hermitePolynomials[i - 2]->getCoefficient(j);
            }

            // compute last two coefficients
            for(int j = numCoeffs - 2; j < numCoeffs; ++j)
            {
                coeffs[j] = hermitePolynomials[i-1]->getCoefficient(j-1);
            }

            // set coefficients of the i-th order Hermite polynomial
            hermitePolynomials[i] = Poly(coeffs);
        }


        // Compute norms of the polynomials
        for(int i = 0; i < maxOrder + 1; ++ i)
        {
            hermiteSquaredNorm[i] = (typeRNum) factorial(i);
        }
    }


    PolynomialConstPtr HermitePolynomialGenerator::getPolynomial(typeInt order) const
    {
        return hermitePolynomials[order];
    }


    typeRNum HermitePolynomialGenerator::getSquaredNorm(typeInt order) const
    {
        return hermiteSquaredNorm[order];
    }

    typeInt HermitePolynomialGenerator::getMaximumOrder() const
    {
        return (typeInt) hermitePolynomials.size() - 1;
    }
}