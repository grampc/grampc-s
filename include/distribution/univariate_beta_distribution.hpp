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


#ifndef UNIVARIATE_BETA_DISTRIBUTION_HPP
#define UNIVARIATE_BETA_DISTRIBUTION_HPP

#include <random>
#include "distribution.hpp"

namespace grampc
{
    // Univariate beta distribution p(x) = 1 / B(p, q) * x^(p-1) * (1 - x)^(q - 1) with 0 < x < 1, p > 0, q > 0
    class UnivariateBetaDistribution : public Distribution
    {
    public:
        // Constructor for the univariate distribution 
        UnivariateBetaDistribution(typeRNum p, typeRNum q);

        // Sample a random numbers of the distribution using the specified random number generator
        virtual const Vector& sample(RandomNumberGenerator& rng) const override;
    
    private:
        mutable typeRNum temp_;
        mutable std::gamma_distribution<typeRNum> distGamma1_;
        mutable std::gamma_distribution<typeRNum> distGamma2_;
        mutable Vector sample_;
    };
    
    // Constructor for the univariate distribution 
    DistributionPtr Beta(typeRNum p, typeRNum q);
}

#endif // UNIVARIATE_BETA_DISTRIBUTION_HPP