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


#ifndef SYMMETRIC_CONSTRAINT_APPROXIMATION
#define SYMMETRIC_CONSTRAINT_APPROXIMATION

#include "chance_constraint_approximation.hpp"

namespace grampc
{
    // Chance constraint approximation assuming a symmetric distribution of the constraints
    class SymmetricConstraintApproximation : public ChanceConstraintApproximation
    {
    public:
        // Constructor for a constraint vector
        SymmetricConstraintApproximation(const Vector& probabilities);

        // Set vector of chance constraint probabilities
        virtual void setConstraintProbability(const Vector& probabilities) override;

        // Get vector of chance constraint probabilities
        virtual const Vector& constraintProbability() const override;

        // Get coefficient z for the constraint tightening of the form z * sqrt(Var{h}) + E{h}
        virtual const Vector& tighteningCoefficient() const override;

    private:
        Vector probabilities_;
        Vector coefficients_;
    };

    // Constructor for a constraint vector
    ChanceConstraintApproximationPtr Symmetric(const Vector& probabilities);
}


#endif // SYMMETRIC_CONSTRAINT_APPROXIMATION