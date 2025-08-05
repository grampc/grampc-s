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


#ifndef CHANCE_CONSTRAINT_APPROXIMATION
#define CHANCE_CONSTRAINT_APPROXIMATION

#include "grampc_interface/grampc_interface.hpp"
#include "util/grampc_s_constants.hpp"

namespace grampc
{
    // Approximation of a chance constraint of the form z * sqrt(Var{h}) + E{h}
    class ChanceConstraintApproximation
    {
    public:
        virtual ~ChanceConstraintApproximation() {};

        // Set vector of chance constraint probabilities
        virtual void setConstraintProbability(const Vector& probabilities) = 0;

        // Get vector of chance constraint probabilities
        virtual const Vector& constraintProbability() const = 0;

        // Get coefficient z for the constraint tightening of the form z * sqrt(Var{h}) + E{h}
        virtual const Vector& tighteningCoefficient() const = 0;
    };

    // Alias
    typedef std::shared_ptr<ChanceConstraintApproximation> ChanceConstraintApproximationPtr;
    typedef std::shared_ptr<const ChanceConstraintApproximation> ChanceConstraintApproximationConstPtr;
}

#endif // CHANCE_CONSTRAINT_APPROXIMATION