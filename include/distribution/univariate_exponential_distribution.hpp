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


#ifndef UNIVARIATE_EXPONENTIAL_DISTRIBUTION_HPP
#define UNIVARIATE_EXPONENTIAL_DISTRIBUTION_HPP

#include <random>
#include "distribution.hpp"

namespace grampc
{
    // Univariate exponential distribution p(x) = lambda * exp(-lambda * x) with lambda > 0, x >= 0
    class UnivariateExponentialDistribution : public Distribution
    {
    public:
        // Constructor for the univariate distribution 
        UnivariateExponentialDistribution(typeRNum lambda);

        // Sample a random numbers of the distribution using the specified random number generator
        virtual const Vector& sample(RandomNumberGenerator& rng) const override;
    
    private:
        mutable std::exponential_distribution<typeRNum> dist_;
        mutable Vector sample_;
    };
    
    // Constructor for the univariate distribution 
    DistributionPtr Exponential(typeRNum lambda);
}

#endif // UNIVARIATE_EXPONENTIAL_DISTRIBUTION_HPP