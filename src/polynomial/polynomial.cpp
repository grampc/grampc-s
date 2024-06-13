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
#include "polynomial/polynomial.hpp"

namespace grampc
{
    Polynomial::Polynomial()
        : numCoefficients_(0)
    {
    }

    Polynomial::Polynomial(const Vector& coefficients)
        : numCoefficients_(coefficients.size()),
          coefficients_(coefficients)
    {
    }

    typeRNum Polynomial::evaluate(ctypeRNum evaluationPoint) const
    {
        return Eigen::poly_eval(coefficients_, evaluationPoint);
    }

    typeRNum Polynomial::gradient(ctypeRNum evaluationPoint) const
    {
        typeRNum out;
        if (numCoefficients_ < 2)
        {
            out = 0.0;
        }
        else
        {
            out = coefficients_[1];
            typeRNum powEval = 1.0;
            for (typeInt i = 2; i < numCoefficients_; ++i)
            {
                powEval *= evaluationPoint;
                out += i * coefficients_[i] * powEval;
            }
        }
        return out;
    }

    typeRNum Polynomial::hessian(ctypeRNum evaluationPoint) const
    {
        if (numCoefficients_ < 3)
        {
            return 0.0;
        }
        else
        {
            typeRNum out = 2.0 * coefficients_[2];
            typeRNum powEval = 1.0;
            for (typeInt i = 3; i < numCoefficients_; ++i)
            {
                powEval *= evaluationPoint;
                out += i * (i - 1) * coefficients_[i] * powEval;
            }
            return out;
        }
    }

    void Polynomial::addPolynomial(const Polynomial& polynomial)
    {
        // Coefficients of the second polynomial
        const Vector& coeff = polynomial.getCoefficients();

        // Polynomial order difference
        typeInt numCoeffDiff = polynomial.numCoefficients_ - numCoefficients_;

        // Increase the number of coefficients if required
        if (numCoeffDiff > 0)
        {
            typeInt elementsOld = coefficients_.rows();
            coefficients_.conservativeResize(elementsOld + numCoeffDiff);
            coefficients_.segment(elementsOld, numCoeffDiff) = Vector::Zero(numCoeffDiff);
            numCoefficients_ = coefficients_.size();
        }

        // Add the coefficients
        coefficients_ += coeff;
    }

    void Polynomial::subtractPolynomial(const Polynomial& polynomial)
    {
        // Coefficients of the second polynomial
        const Vector& coeff = polynomial.getCoefficients();

        // Polynomial order difference
        typeInt numCoeffDiff = polynomial.numCoefficients_ - numCoefficients_;

        // Increase the number of coefficients if required
        if (numCoeffDiff > 0)
        {
            typeInt elementsOld = coefficients_.rows();
            coefficients_.conservativeResize(elementsOld + numCoeffDiff);
            coefficients_.segment(elementsOld, numCoeffDiff) = Vector::Zero(numCoeffDiff);
            numCoefficients_ = coefficients_.size();
        }

        // Add the coefficients
        coefficients_ -= coeff;
    }

    void Polynomial::multiplyPolynomial(const Polynomial& polynomial)
    {
        const Vector& coeff_2 = polynomial.getCoefficients();
        typeInt numCoeff_2 = polynomial.getNumCoefficients();

        // multiply each coefficient with all coefficients of the second polynomial
        Vector tempCoef = Vector::Zero(numCoefficients_ + numCoeff_2 - 1);
        for (typeInt i = 0; i < numCoefficients_; ++i)
        {
            for (typeInt j = 0; j < numCoeff_2; ++j)
            {
                tempCoef[i + j] += coefficients_[i] * coeff_2[j];
            }
        }

        // Set Coefficients
        coefficients_ = tempCoef;
        numCoefficients_ = coefficients_.size();
    }

    void Polynomial::multiplyScalar(const typeRNum factor)
    {
        coefficients_ *= factor;
    }

    void Polynomial::setCoefficients(const Vector& coefficients)
    {
        coefficients_ = coefficients;
        numCoefficients_ = coefficients_.size();
    }

    const Vector& Polynomial::getCoefficients() const
    {
        return coefficients_;
    }

    typeRNum Polynomial::getCoefficient(int index) const
    {
        return coefficients_[index];
    }

    typeInt Polynomial::getNumCoefficients() const
    {
        return numCoefficients_;
    }

    PolynomialPtr Poly()
    {
        return PolynomialPtr(new Polynomial());
    }

    PolynomialPtr Poly(const Vector& coefficients)
    {
        return PolynomialPtr(new Polynomial(coefficients));
    }

}