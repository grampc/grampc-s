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


#include "constraint_approx/chebyshev_constraint_approximation.hpp"

namespace grampc
{
    ChebyshevConstraintApproximation::ChebyshevConstraintApproximation(const Vector& probabilities)
    : probabilities_(probabilities),
      coefficients_(probabilities.rows())
    {
        // Go through all constraints
        for(typeInt i = 0; i < probabilities_.size(); ++i)
        {
            // Set coefficient
            coefficients_[i] = std::sqrt(probabilities_[i] / (1.0 - probabilities_[i]));
        }
    }


    void ChebyshevConstraintApproximation::setConstraintProbability(const Vector& probabilities)
    {
        probabilities_ = probabilities;

        // Number of constraints
        typeInt numConstraints = probabilities_.rows();

        // Resize coefficients
        coefficients_.resize(numConstraints);

        // Go through all constraints
        for(typeInt i = 0; i < numConstraints; ++i)
        {
            // Set coefficient
            coefficients_[i] = std::sqrt(probabilities_[i] / (1.0 - probabilities_[i]));
        }

    }


    const Vector& ChebyshevConstraintApproximation::constraintProbability() const
    {
        return probabilities_;
    }


    const Vector& ChebyshevConstraintApproximation::tighteningCoefficient() const
    {
        return coefficients_;
    }

    ChanceConstraintApproximationPtr Chebyshev(const Vector& probabilities)
    {
        return ChanceConstraintApproximationPtr(new ChebyshevConstraintApproximation(probabilities));
    }
}