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


#include "polynomial/multivariate_polynomial.hpp"

namespace grampc
{
    MultivariatePolynomial::MultivariatePolynomial(const std::vector<PolynomialConstPtr>& univariatePolynomials, const std::vector<typeRNum>& univariateSquaredNorm)
        : polynomials_(univariatePolynomials),
          numVariables_(univariatePolynomials.size()),
          squaredNorm_(1.0)
    {
        // Compute squared norm
        for (typeInt i = 0; i < numVariables_; ++i)
        {
            squaredNorm_ *= univariateSquaredNorm[i];
        }
    }

    typeRNum MultivariatePolynomial::evaluate(const Vector& evaluationPoint) const
    {
        typeRNum out = 1.0;

        for (typeInt i = 0; i < numVariables_; ++i)
        {
            out *= polynomials_[i]->evaluate(evaluationPoint[i]);
        }
        return out;
    }

    void MultivariatePolynomial::setPolynomials(std::vector<PolynomialConstPtr>& univariatePolynomials, std::vector<typeRNum>& univariateSquaredNorm)
    {
        polynomials_ = univariatePolynomials;
        numVariables_ = polynomials_.size();

        // Compute squared norm
        squaredNorm_ = 1.0;
        for (typeInt i = 0; i < numVariables_; ++i)
        {
            squaredNorm_ *= univariateSquaredNorm[i];
        }
    }

    const std::vector<PolynomialConstPtr> &MultivariatePolynomial::polynomials() const
    {
        return polynomials_;
    }

    PolynomialConstPtr MultivariatePolynomial::polynomials(int index) const
    {
        return polynomials_[index];
    }

    typeInt MultivariatePolynomial::numVariables() const
    {
        return numVariables_;
    }

    typeRNum MultivariatePolynomial::squaredNorm() const
    {
        return squaredNorm_;
    }

    MultivariatePolynomialPtr MultivarPoly(const std::vector<PolynomialConstPtr>& univariatePolynomials, const std::vector<typeRNum>& univariateSquaredNorm)
    {
        return MultivariatePolynomialPtr(new MultivariatePolynomial(univariatePolynomials, univariateSquaredNorm));
    }
}