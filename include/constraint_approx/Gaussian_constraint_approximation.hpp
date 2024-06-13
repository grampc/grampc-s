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


#ifndef GAUSSIAN_CONSTRAINT_APPROXIMATION
#define GAUSSIAN_CONSTRAINT_APPROXIMATION

#include "chance_constraint_approximation.hpp"

namespace grampc
{
    // Chance constraint approximation assuming a Gaussian distribution of the constraints
    class GaussianConstraintApproximation : public ChanceConstraintApproximation
    {
    public:
        // Constructor for a constraint vector
        GaussianConstraintApproximation(const Vector& probabilities);

        // Set vector of chance constraint probabilities
        virtual void setConstraintProbability(const Vector& probabilities) override;

        // Get vector of chance constraint probabilities
        virtual const Vector& constraintProbability() const override;

        // Get coefficient z for the constraint tightening of the form z * sqrt(Var{h}) + E{h}
        virtual const Vector& tighteningCoefficient() const override;

    private:
        Vector probabilities_;
        Vector coefficients_;

        const Eigen::Matrix<typeRNum, 200, 1> integral_vec_;
        const Eigen::Matrix<typeRNum, 200, 1> quantile_vec_;
    };

    // Constructor for a constraint vector
    ChanceConstraintApproximationPtr GaussianApprox(const Vector& probabilities);
}


#endif // GAUSSIAN_CONSTRAINT_APPROXIMATION